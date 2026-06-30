###############################################################################
## 01_generate_ipw_weights.R
##
## Generates the loss-to-follow-up inverse-probability weights used throughout
## the secondary HTE analysis. For each of the 13 arms it fits a SuperLearner
## ensemble predicting 4-week follow-up completion from baseline covariates,
## then forms stabilised inverse-probability-of-completion weights:
##
##     raw_w  = 1 / p_hat   (completers only; 0 otherwise)
##     w_norm = raw_w * n_completers / sum(raw_w)
##
## A cross-validated SuperLearner is also fit per arm to record the CV-AUC of
## the completion model (reported in output/tables/loss_to_follow_up.csv).
##
## Outputs:
##   data/sl_weights.csv          one row per participant (arm, pid, w_norm)
##   output/tables/loss_to_follow_up.csv   per-arm CV-AUC + weight percentiles
##
## Run standalone:  Rscript R/01_generate_ipw_weights.R
## Or source() from the main analysis document.
##
## NOTE: This step is computationally intensive (a full + cross-validated
## SuperLearner ensemble per arm). Expect a long run time.
###############################################################################

suppressMessages({
  library(SuperLearner)
  library(glmnet)
  library(ranger)
  library(cvAUC)
  library(caret)
  library(dplyr)
  library(tibble)
})

## Resolve paths relative to this script's project root (the directory that
## contains data/ and output/). Works whether sourced or run via Rscript.
if (!exists("REPRO_ROOT")) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  REPRO_ROOT <- if (length(file_arg) == 1) {
    normalizePath(file.path(dirname(file_arg), ".."))
  } else {
    normalizePath(".")
  }
}
source(file.path(REPRO_ROOT, "R", "sl_helpers.R"))

## ---- Data ------------------------------------------------------------------
mydata <- read.csv(file.path(REPRO_ROOT, "data", "ml_clean.csv"), header = TRUE)
mydata <- mydata[, -1]                                   # drop row-number column
mydata <- mydata %>% select(-help_prefer_not_say_baseline)  # zero-variance
mydata <- mydata %>% mutate(fup_complete = ifelse(!is.na(phq9_score_followup), 1, 0))

## ---- SuperLearner configuration --------------------------------------------
simple_SL_library <- build_ipw_sl_library()
cvControl <- SuperLearner.CV.control(V = 5L, stratifyCV = TRUE)

get_cv_auc <- function(y, cv_pred, valid_rows) {
  if (length(unique(y)) < 2) return(NA_real_)
  out <- try(cvAUC::ci.cvAUC(cv_pred, y, folds = valid_rows), silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else unname(out$cvAUC)
}

probs  <- c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1)
labels <- c("Min", "1st percentile", "5th percentile", "10th percentile",
            "25th percentile", "50th percentile", "75th percentile",
            "90th percentile", "95th percentile", "99th percentile", "Max")

arms <- sort(unique(mydata$group_numeric))
summary_rows <- list()
weights_long <- list()

for (arm in arms) {
  message(sprintf("[IPW] arm %d", arm))
  dat_arm <- mydata %>% filter(group_numeric == arm)

  n_total     <- nrow(dat_arm)
  n_missing   <- sum(is.na(dat_arm$phq9_score_followup))
  pct_missing <- if (n_total > 0) n_missing / n_total else NA_real_

  X <- as.data.frame(dat_arm[, 4:102])                 # baseline covariates
  y <- ifelse(is.na(dat_arm$phq9_score_followup), 0, 1)

  ## Cross-validated SuperLearner -> CV-AUC of the completion model
  set.seed(12345)
  cv_fit <- CV.SuperLearner(
    Y = y, X = X, family = binomial(),
    SL.library = simple_SL_library, verbose = FALSE,
    cvControl = cvControl, innerCvControl = list(cvControl)
  )
  cv_auc <- get_cv_auc(y, cv_fit$SL.predict, cv_fit$validRows)

  ## Full SuperLearner -> predicted completion probability p_hat
  set.seed(12345)
  sl_fit <- SuperLearner(
    Y = y, X = X, family = binomial(),
    SL.library = simple_SL_library, verbose = FALSE, cvControl = cvControl
  )
  pred <- as.numeric(sl_fit$SL.predict)

  ## Stabilised inverse-probability-of-completion weights
  raw_w  <- ifelse(y == 1, 1 / pred, 0)
  denom  <- sum(raw_w[y == 1])
  w_norm <- ifelse(y == 1 & denom > 0, raw_w * sum(y == 1) / denom, 0)

  weights_long[[as.character(arm)]] <- tibble(
    intervention = arm, pid = dat_arm$pid,
    fup_complete = y,   w_norm = w_norm
  )

  w_comp <- w_norm[y == 1]
  qs <- if (length(w_comp) > 0)
    as.numeric(quantile(w_comp, probs = probs, na.rm = TRUE, type = 7))
  else rep(NA_real_, length(probs))

  arm_row <- tibble(
    intervention = arm, n_missing = n_missing,
    pct_missing  = round(pct_missing, 4),
    auc          = ifelse(is.na(cv_auc), NA_real_, round(cv_auc, 4))
  )
  names(qs) <- labels
  summary_rows[[as.character(arm)]] <- bind_cols(arm_row, as_tibble_row(round(qs, 4)))
}

summary_table   <- bind_rows(summary_rows)
weights_long_df <- bind_rows(weights_long)

dir.create(file.path(REPRO_ROOT, "output", "tables"),
           showWarnings = FALSE, recursive = TRUE)
write.csv(weights_long_df, file.path(REPRO_ROOT, "data", "sl_weights.csv"))
write.csv(summary_table,
          file.path(REPRO_ROOT, "output", "tables", "loss_to_follow_up.csv"))

message("[IPW] done. Wrote data/sl_weights.csv and output/tables/loss_to_follow_up.csv")
print(summary_table)
