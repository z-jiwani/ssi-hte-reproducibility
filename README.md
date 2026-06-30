# Reproducibility package — Heterogeneous treatment effects of single-session interventions for depression

This package reproduces, from the cleaned analytic data and **without any
cached intermediates**, the secondary heterogeneous-treatment-effect (HTE)
analysis of a crowdsourced megastudy of digital single-session interventions
(SSIs) for depression (N = 7,505 randomized to 12 SSIs or a passive control;
4-week PHQ-9 follow-up).

The analysis (1) estimates individual-level conditional average treatment
effects (CATEs) for each SSI versus control with generalized causal forests,
(2) clusters participants on their 12-dimensional CATE profiles with latent
profile analysis, (3) tests cluster-specific intervention effects, and
(4) compares assignment strategies by bootstrap.

## Contents

```
.
├── README.md
├── ssi_hte_analysis.qmd          # main analysis (renders to HTML)
├── ssi_hte_analysis.html         # rendered analysis with all results (provided)
├── R/
│   ├── sl_helpers.R              # serpentine folds + SuperLearner screeners/learners
│   └── 01_generate_ipw_weights.R # loss-to-follow-up IPW weight generation
├── data/
│   ├── ml_clean.csv              # cleaned analytic dataset (N = 7,505)
│   ├── ssi_characteristics.csv   # SSI design-feature ratings
│   └── sl_weights.csv            # loss-to-follow-up weights (produced by step 1)
└── output/
    ├── tables/                   # all tables written here
    └── figures/                  # all figures written here
```

## Data

The cleaned analytic dataset (`data/ml_clean.csv`) is derived from the publicly
available crowdsourced megastudy of single-session interventions for depression
(de-identified; <https://osf.io/agvh6>). Each row is a participant; columns are
a participant id, the randomized arm (`group_numeric`, 1–12 = SSIs,
13 = passive control), the 4-week PHQ-9 outcome, and baseline covariates.
`data/ssi_characteristics.csv` holds expert design-feature ratings for the 12
SSIs. Only the variables used by this analysis are included.

## Software

- R (developed under R 4.4.1)
- Quarto (to render the analysis document)
- R packages: `tidyverse`, `grf`, `mclust`, `glmnet`, `tidymodels`, `emmeans`,
  `broom`, `boot`, `corrplot`, `reshape2`, `ggrepel`, `car` (analysis) and
  `SuperLearner`, `ranger`, `cvAUC`, `caret` (IPW weights).

The serpentine cross-validation-fold routine is included directly in
`R/sl_helpers.R`, so no non-CRAN package is required.

## How to reproduce

Run from this directory.

**Step 1 — loss-to-follow-up weights** (computationally intensive; fits a
cross-validated SuperLearner completion model per arm):

```bash
Rscript R/01_generate_ipw_weights.R
```

This writes `data/sl_weights.csv` and `output/tables/loss_to_follow_up.csv`.

**Step 2 — full analysis** (regenerates the causal forests, the latent profile
analysis, the LASSO characterization models, and the 5,000-resample bootstrap
from scratch):

```bash
quarto render ssi_hte_analysis.qmd
```

This produces `ssi_hte_analysis.html` and writes every table to
`output/tables/` and every figure to `output/figures/`.

> A pre-rendered `ssi_hte_analysis.html` is included so results can be inspected
> without re-running the pipeline.

## Reproducibility notes

- All stochastic steps are seeded, so re-runs reproduce the reported values
  (the IPW weights, CATE matrix, LPA solution, LASSO models, and bootstrap
  contrasts) up to floating-point tolerance on the same platform.
- The 4-week PHQ-9 outcome is reverse-scored at load (`27 - PHQ-9`) so that
  higher values denote improvement; treatment-effect coefficients are therefore
  positive when the intervention reduces depressive symptoms.
- Loss to follow-up is handled by stabilized inverse-probability-of-completion
  weights applied throughout CATE estimation and all outcome models.

## Data source

The analytic dataset (`data/ml_clean.csv`) is derived from the publicly
available, de-identified trial dataset hosted on the Open Science Framework
(<https://osf.io/agvh6>). Only the variables used by this analysis are
included here.
