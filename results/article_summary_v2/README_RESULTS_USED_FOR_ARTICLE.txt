These are the numerical results selected for the revised manuscript.

Use this folder for article tables:
  results/article_summary_v2

Use for estimation/RMSE/sample-size summaries:
  summary_estimation_article.csv
  summary_fit_quality_article.csv

Use for LR tests:
  summary_lr_lambda_article.csv
  summary_lr_delta_article.csv

Use only as diagnostic:
  summary_wald_coverage_diagnostic.csv

Do not use for article estimation tables:
  results/summary/summary_estimation.csv
  results/article_summary_current/summary_estimation_article.csv

Reason:
  the earlier estimation summaries were produced with weaker optimization or with an overly restrictive strict-success filter. The v2 summaries use the reinforced optimization run and all valid maximum-likelihood fits.

Inference note:
  Wald coverage for skewness components is diagnostic only. Formal inference in the manuscript should rely on likelihood-ratio tests, with parametric-bootstrap calibration for the boundary test on the warp exponents.
