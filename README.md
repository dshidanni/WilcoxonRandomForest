# Wilcoxon Random Forests for Robust Distributional Prediction

*Author: Danni Shi, danni.shi[at]vanderbilt.edu; dshidanni[at]outlook.com*

This repository contains all the R codes accompanying the manuscript **"Wilcoxon Random Forests for Robust Distributional Prediction"**.

We propose a Wilcoxon random forest (WRF) for prediction of outcome distributions. The method grows trees using a rank-based splitting criterion that targets stochastic separation between the outcome distributions in child nodes and estimates conditional distributions using forest-weighted empirical CDFs. The WRF can be used to estimate conditional quantiles, means, and threshold probabilities. Because the splitting criterion depends only on outcome ranks, WRF is invariant to strictly increasing transformations and can accommodate continuous, ordinal, and mixed outcomes, including biomarkers with detection limits. We establish consistency, we investigate performance with simulations, and we illustrate the method by predicting HIV outcomes in a multicenter Latin American cohort.

This repository is *not an R package* as of September 25, 2026, and the public interface may change. A package is planned after manuscript submission.

## `R/`

Source implementation of:

-   Wilcoxon, SSE and pinball-loss regression tree splitting;
-   WRF and quantile regression forest fitting;
-   conditional quantiles, means, medians, threshold probabilities, ECDFs, and forest weights;
-   CRPS, Cramér–von Mises calibration;
-   out-of-bag prediction and CRPS-based permutation variable importance;
-   optional CPM refinement of fitted Wilcoxon trees/forests.

See [`R/README.md`](R/README.md) for the API, estimator conventions, examples, and limitations.

## `Simulations/`

Scripts for Section 4 and Supplementary Section S2, including:

-   a seven-region tree-structured data-generating process;
-   its strictly monotone right-skewing transformation;
-   a smooth nonlinear model with right-skewed heteroscedastic errors;
-   single-tree and forest comparisons;
-   creation of manuscript and supplementary result tables.

See [`Simulations/README.md`](Simulations/README.md) for the full design, run order, expected outputs, and computing requirements.

## `Application/`

Public workflow for the HIV application using the synthetic data. The confidential clinical data are not included. The synthetic data and scripts demonstrate:

-   6-month CD4 distribution prediction and threshold probabilities `P(CD4 <= 200, 350, 500 | X)`;
-   6-month viral suppression probability, `P(VL < 80 | X)`;
-   comparison with logistic regression and classification random forest;
-   cross-validated Brier score, CRPS, calibration, and predicted-CDF displays.

The synthetic results are illustrative and will not numerically reproduce the confidential-data results in the manuscript, but the R scripts the R codes execute the exact same pipeline in Section 5 and Supplementary S3 of the manuscript. See [`Application/README.md`](Application/README.md).


## Computation details

All code was implemented in R 4.5.3. We also use CRAN packages

-   `rms` package (version 8.1-0 for CPM-based approaches) for CPM refinement. For CPM refinement, please make sure that the `rms` package version is at least 7.0-0.
-   `Matrix` package (version 1.7-4) for matrix computations.
-   `synthpop` package (version 1.9-2) for generating synthetic datasets.

## Reference

Shi D, Shepherd BE, Li C, (2026) Wilcoxon Random Forests for Robust Distributional Prediction.

Harrell Jr FE (2025). _rms: Regression Modeling Strategies_. doi:10.32614/CRAN.package.rms <https://doi.org/10.32614/CRAN.package.rms>, R package version 8.1-0, <https://CRAN.R-project.org/package=rms>.

Bates D, Maechler M, Jagan M (2025). _Matrix: Sparse and Dense Matrix Classes and Methods_. doi:10.32614/CRAN.package.Matrix <https://doi.org/10.32614/CRAN.package.Matrix>, R package version 1.7-4, <https://CRAN.R-project.org/package=Matrix>.

Beata Nowok, Gillian M. Raab, Chris Dibben (2016). synthpop: Bespoke Creation of Synthetic Data in R. Journal of Statistical Software, 74(11), 1-26. doi:10.18637/jss.v074.i11
