###############################################################################
## sl_helpers.R
##
## Helper functions for the secondary HTE analysis of the crowdsourced
## megastudy of single-session interventions (SSIs) for depression.
##
## Contents:
##   1. Serpentine sort (self-contained; used for cross-validation folds)
##   2. SuperLearner screeners and custom learners (used for the loss-to-
##      follow-up inverse-probability weights)
##   3. create.Learner.grid() and SL.ranger.1() learner constructors
##
## This file is sourced by 01_generate_ipw_weights.R and by the main
## analysis document. It has no external file-path dependencies.
###############################################################################

## ===========================================================================
## 1. Serpentine sort
## ---------------------------------------------------------------------------
## Stratified, balanced ordering used to build cross-validation folds. The
## data are repeatedly split on the leading variable and sorted in
## alternating (up / down) order on the next variable, so that adjacent rows
## are balanced across the sort variables. Folds are then assigned by taking
## every k-th row of the serpentine-sorted data.
##
## (Self-contained reimplementation; produces row orderings identical to the
## `samplingTools::serpentine` function used during the original analysis.)
## ===========================================================================

single_serp <- function(data, var_desc, var_serp) {
  split_data <- split(data, dplyr::pull(data, !!var_desc))
  split_data[c(TRUE, FALSE)] <- lapply(split_data[c(TRUE, FALSE)], function(x) {
    dplyr::arrange(x, !!var_serp)
  })
  split_data[c(FALSE, TRUE)] <- lapply(split_data[c(FALSE, TRUE)], function(x) {
    dplyr::arrange(x, dplyr::desc(!!var_serp))
  })
  dplyr::bind_rows(split_data)
}

serpentine <- function(data = NULL, ...) {
  orig_names <- names(data)
  full_vars  <- rlang::quos(...)
  if (length(full_vars) < 2) stop("Need at least two variables to perform sort")

  data <- single_serp(data, full_vars[[1]], full_vars[[2]])
  data <- dplyr::mutate(data, row_n = dplyr::row_number())
  data <- dplyr::group_by(data, !!!full_vars[1:2])
  data <- dplyr::mutate(data, new_group_num = min(row_n))
  data <- dplyr::ungroup(data)

  if (length(full_vars) >= 3) {
    for (j in 3:length(full_vars)) {
      data <- single_serp(data, rlang::quo(new_group_num), full_vars[[j]])
      data <- dplyr::mutate(data, row_n = dplyr::row_number())
      data <- dplyr::group_by(data, !!!full_vars[1:j])
      data <- dplyr::mutate(data, new_group_num = min(row_n))
      data <- dplyr::ungroup(data)
    }
  }
  dplyr::select(data, dplyr::one_of(orig_names))
}

## Serpentine-sort within a stratum, then assign `nfolds` folds.
## Sort vars: outcome, baseline PHQ-9, follow-up indicator, IPW weight.
make_folds_serpentine <- function(dat, nfolds = 10, seed = 123) {
  set.seed(seed)
  start <- sample(0:(nfolds - 1), 1)
  dat_sorted <- do.call(
    serpentine,
    c(list(data = dat),
      lapply(c("phq9_score_followup", "phq9_pre_raw", "fup_complete", "w_norm"),
             rlang::sym))
  )
  dat_sorted$fold <- (((seq_len(nrow(dat_sorted)) - 1 + start) %% nfolds) + 1L)
  dat_sorted
}


## ===========================================================================
## 2. SuperLearner learner constructors
## ===========================================================================

## create.Learner.grid(): like SuperLearner::create.Learner() but accepts a
## tuning grid in which each row sets a full combination of hyper-parameters.
create.Learner.grid <- function(base_learner, params = list(), tune = list(),
                                 tunegrid = data.frame(),
                                 env = parent.frame(), name_prefix = base_learner,
                                 detailed_names = FALSE, verbose = FALSE) {
  if (length(tunegrid) > 0) {
    tuneGrid <- tunegrid; names <- rep("", nrow(tuneGrid)); max_runs <- nrow(tuneGrid)
  } else if (length(tune) > 0) {
    tuneGrid <- expand.grid(tune, stringsAsFactors = FALSE)
    names <- rep("", nrow(tuneGrid)); max_runs <- nrow(tuneGrid)
  } else {
    max_runs <- 1; tuneGrid <- NULL; names <- c()
  }
  for (i in seq(max_runs)) {
    name <- paste(name_prefix, i, sep = "_")
    if (length(tuneGrid) > 0) {
      g <- tuneGrid[i, , drop = FALSE]
      if (detailed_names) name <- do.call(paste, c(list(name_prefix), g, list(sep = "_")))
    } else g <- c()
    names[i] <- name
    fn_params <- ""
    all_params <- c(as.list(g), params)
    for (name_i in names(all_params)) {
      val <- all_params[[name_i]]
      if (!is.null(val) && val != "NULL") {
        if (class(val) == "character") val <- paste0('"', val, '"')
        fn_params <- paste0(fn_params, ", ", name_i, "=", val)
      }
    }
    fn <- paste0(name, " <- function(...) ", base_learner, "(...", fn_params, ")")
    if (verbose) cat(fn, "\n")
    eval(parse(text = fn), envir = env)
  }
  invisible(list(grid = tuneGrid, names = names,
                 base_learner = base_learner, params = params))
}

## SL.ranger.1(): SL.ranger wrapper that exposes max.depth.
SL.ranger.1 <- function(Y, X, newX, family, obsWeights,
                        num.trees = 500, mtry = floor(sqrt(ncol(X))),
                        write.forest = TRUE,
                        probability = family$family == "binomial",
                        min.node.size = ifelse(family$family == "gaussian", 5, 1),
                        max.depth = NULL, replace = TRUE,
                        sample.fraction = ifelse(replace, 1, 0.632),
                        num.threads = 1, verbose = TRUE, ...) {
  SuperLearner:::.SL.require("ranger")
  if (family$family == "binomial") Y <- as.factor(Y)
  if (is.matrix(X)) X <- data.frame(X)
  fit <- ranger::ranger(`_Y` ~ ., data = cbind("_Y" = Y, X),
                        num.trees = num.trees, mtry = mtry,
                        min.node.size = min.node.size, max.depth = max.depth,
                        replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest,
                        probability = probability, num.threads = num.threads,
                        verbose = verbose)
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") pred <- pred[, "1"]
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  list(pred = pred, fit = fit)
}


## ===========================================================================
## 3. SuperLearner screeners (two-part: rare -> correlation -> LASSO / ranger)
## ===========================================================================

screen.lasso.sub <- function(Y, X, family, subs = c(1:ncol(X)), alpha = 1,
                             minscreen = 10, nfolds = 10, nlambda = 100,
                             obsWeights, dfmax = ncol(X) + 1, ...) {
  SuperLearner:::.SL.require("glmnet")
  if (!is.matrix(X)) X <- model.matrix(~ -1 + ., X)
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "default",
                             nfolds = nfolds, weights = obsWeights, dfmax = dfmax,
                             family = family$family, alpha = alpha, nlambda = nlambda, ...)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  if (sum(whichVariable) < minscreen) {
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
  }
  whichVariable
}

screen.ranger.sub <- function(Y, X, family, subs = c(1:ncol(X)), nVar = ncol(X),
                              ntree = 1000, mtry = 34, maxdepth = 5,
                              splitrule = ifelse(family$family == "gaussian", "variance", "gini"),
                              obsWeights, ...) {
  X_ranger <- X
  SuperLearner:::.SL.require("ranger")
  if (family$family == "gaussian") {
    rank.rf.fit <- ranger::ranger(Y ~ ., data = X_ranger, num.trees = ntree,
                                  splitrule = splitrule, max.depth = maxdepth,
                                  mtry = ifelse(mtry < ncol(X), mtry, ),
                                  importance = "impurity_corrected",
                                  case.weights = obsWeights)
  }
  if (family$family == "binomial") {
    rank.rf.fit <- ranger::ranger(as.factor(Y) ~ ., data = X_ranger, num.trees = ntree,
                                  splitrule = splitrule, max.depth = maxdepth,
                                  mtry = ifelse(mtry < ncol(X), mtry, ),
                                  importance = "impurity_corrected",
                                  case.weights = obsWeights)
  }
  as.vector(rank(-rank.rf.fit$variable.importance[1:ncol(X_ranger)]) <= nVar)
}

screen.rare <- function(Y, X, family, thresh = 5, ...) {
  if (family$family == "binomial") {
    max_limit <- min(table(Y)) - thresh
    min_Y_val <- names(table(Y))[table(Y) == min(table(Y))]
    X_rare <- X[Y == min_Y_val, ]
  } else if (family$family == "gaussian") {
    max_limit <- nrow(X) - thresh; X_rare <- X
  }
  max_freq <- apply(X_rare, 2, function(x) max(table(x)))
  max_freq <= max_limit
}

screen.corr <- function(Y, X, obsWeights, r = 0.95, ...) {
  corr <- cov.wt(X, wt = obsWeights, cor = TRUE)$cor
  if (any(corr[upper.tri(corr, diag = FALSE)] >= r)) {
    high_corr <- caret::findCorrelation(corr, cutoff = r, names = TRUE)
    !colnames(X) %in% high_corr
  } else rep(TRUE, ncol(X))
}

screen.two_part.lasso <- function(Y, X, thresh = 5, r = 0.95, family, obsWeights,
                                  minscreen = 10, dfmax = ncol(X) + 1, ...) {
  whichNotRare <- screen.rare(Y = Y, X = X, family = family, thresh = thresh, ...)
  X_notRare <- X[, whichNotRare]
  whichCorr <- screen.corr(Y = Y, X = X_notRare, r = r, obsWeights = obsWeights, ...)
  screen.x <- X_notRare[, whichCorr]
  if (length(screen.x) < minscreen) minscreen <- length(screen.x)
  lassoVariables <- screen.lasso.sub(Y = Y, X = screen.x, family = family,
                                     obsWeights = obsWeights, alpha = 1,
                                     minscreen = minscreen, dfmax = dfmax, ...)
  colnames(X) %in% colnames(screen.x)[lassoVariables]
}

screen.two_part.ranger <- function(Y, X, thresh = 5, r = 0.95, family, obsWeights,
                                   nVar = 30, ...) {
  whichNotRare <- screen.rare(Y = Y, X = X, family = family, thresh = thresh, ...)
  X_notRare <- X[, whichNotRare]
  whichCorr <- screen.corr(Y = Y, X = X_notRare, r = r, obsWeights = obsWeights, ...)
  screen.x <- X_notRare[, whichCorr]
  rangerVariables <- screen.ranger.sub(Y = Y, X = screen.x, family = family,
                                       obsWeights = obsWeights, nVar = nVar, ...)
  colnames(X) %in% colnames(screen.x[, rangerVariables])
}

## Screener instances at three covariate-budget levels
nvar1 <- 25; nvar2 <- 40; nvar3 <- 60
screen.two_part.lasso_1  <- function(...) screen.two_part.lasso(dfmax = nvar1, ...)
screen.two_part.lasso_2  <- function(...) screen.two_part.lasso(dfmax = nvar2, ...)
screen.two_part.lasso_3  <- function(...) screen.two_part.lasso(dfmax = nvar3, ...)
screen.two_part.ranger_1 <- function(...) screen.two_part.ranger(nVar = nvar1, ...)
screen.two_part.ranger_2 <- function(...) screen.two_part.ranger(nVar = nvar2, ...)
screen.two_part.ranger_3 <- function(...) screen.two_part.ranger(nVar = nvar3, ...)


## ===========================================================================
## 4. Build the SuperLearner library used for the loss-to-follow-up model
## ===========================================================================

build_ipw_sl_library <- function(env = globalenv()) {
  # NOTE: the individual learner functions (SL.glmnet_0, SL.ranger_*, ...) must
  # be created in an environment SuperLearner can later look them up in (it
  # resolves learner names with get() against the calling environment). We
  # therefore create them in `env` (the global environment by default), not in
  # this function's transient evaluation frame.
  tune.glmnet <- list(alpha = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))

  set.seed(250829)
  n_tunes <- 6; n_real_Y <- 500
  tune.grid.ranger1 <- expand.grid(
    num.trees = 5000L, max.depth = 3:5L,
    min.node.size = as.integer(n_real_Y * c(0.05, 0.1, 0.15)),
    mtry = as.integer(c(max(2, sqrt(nvar1) / 2), sqrt(nvar1), min(nvar1, sqrt(nvar1) * 2)))
  ) |> data.frame()
  selected.tune.ranger1 <- tune.grid.ranger1[sample(nrow(tune.grid.ranger1), n_tunes / 3), ]

  grid.glmnet  <- SuperLearner::create.Learner("SL.glmnet", tune = tune.glmnet,
                                               detailed_names = TRUE, env = env)
  grid.ranger1 <- create.Learner.grid("SL.ranger.1", detailed_names = TRUE,
                                       tunegrid = selected.tune.ranger1,
                                       name_prefix = "SL.ranger", env = env)

  grid.glmnet.screen  <- lapply(grid.glmnet$names, c,
                                "screen.two_part.lasso_1",
                                "screen.two_part.lasso_2",
                                "screen.two_part.lasso_3")
  grid.ranger.screen1 <- lapply(grid.ranger1$names, c, "screen.two_part.ranger_1")

  c(list("SL.mean"), grid.glmnet.screen, grid.ranger.screen1)
}
