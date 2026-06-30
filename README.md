# HTE analysis of single-session interventions for depression

Reproducible code for a secondary heterogeneous-treatment-effect analysis of a
crowdsourced megastudy of digital SSIs for depression (N = 7,505; 12 SSIs vs.
control). Requires R and Quarto.

```bash
Rscript R/01_generate_ipw_weights.R   # step 1: loss-to-follow-up weights
quarto render ssi_hte_analysis.qmd    # step 2: full analysis -> ssi_hte_analysis.html
```

`data/ml_clean.csv` is derived from the de-identified trial dataset on the
Open Science Framework (https://osf.io/agvh6).
