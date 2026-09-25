# R source code: Wilcoxon regression trees and Wilcoxon random forests

This directory includes the implementation of the methods described in the manuscript **“Wilcoxon Random Forests for Robust Distributional Prediction.”** The code fits distributional regression trees and random forests using squared-error, Wilcoxon rank-based, or pinball-loss splitting, and obtains predictions from forest-weighted empirical cumulative distribution functions (CDFs).

> **Note:** This is research code. This is not an R package and has not yet been intended as a stable public API. The authors plan to develop a package after manuscript submission.

## Existing files

| File | Purpose |
|------------------------------------|------------------------------------|
| `myRF.R` | Core implementation. Fits trees and forests with `myRF()`; provides prediction, out-of-bag evaluation, CRPS-based permutation variable importance, and an out-of-bag pinball-loss summary. |
| `CPMtree.R` | Optional cumulative probability model (CPM) refinement of a fitted `myRF` object. Source this file after `myRF.R`. |

The manuscript’s primary forest estimator is the ECDF-based Wilcoxon random forest (WRF). CPM refinement is used for single-tree analyses and an exploratory forest analysis; it is not part of the primary WRF estimator.

All scripts assume that the working directory is the root directory of this repository. File paths are specified relative to the repository root.

## Requirements

The source files require:

-   R;
-   `Matrix`, for sparse forest-weight matrices;
-   `rms` version 7.0-0 or newer, only for `CPMtree.R`.

Install the direct dependencies with:

``` r
install.packages(c("Matrix", "rms"))
```

## Loading the code

From the repository root:

``` r
source("R/myRF.R")

# Only needed for CPM refinement:
source("R/CPMtree.R")
```

`predict.myRF()` and `predict.myCPMForest()` are S3 prediction methods, so either `predict(object, ...)` or the full method name may be used.

The functions `predict_oob.myRF()`, `variable_importance_crps.myRF()`, and `oob_pinball_loss.myRF()` are not S3 generics and should be called by their full names.

## Data requirements

### Response

The response must be represented numerically and its numeric order must have the intended scientific meaning. This includes continuous, count, ordinal, and mixed discrete–continuous outcomes.

The current implementation applies `as.numeric()` to the response internally. Therefore:

-   do not pass an unordered factor or character response;
-   convert an ordered factor explicitly to meaningful numeric scores before fitting;
-   for an outcome subject to a detection limit, use a numeric code that preserves the ordering. For example, a common “below limit” code must be strictly below every measured value when it is intended to represent the lowest ordered category.

### Predictors

Predictors may be numeric, factor, or character variables. Numeric splits have the form `x <= s`. Categorical splits send a stored set of levels to the left child. A categorical level not observed during fitting is sent to the right child by the current prediction code.

### Missing and infinite values

Missing values are not imputed. `myRF()` stops when the model frame contains `NA`, and prediction stops when `newdata` contains `NA`. Imputation or complete-case processing must therefore be completed before fitting.

For SSE and pinball splitting, infinite numeric values cause an error. Wilcoxon splitting can rank infinite values, but the code gives a warning because they usually indicate a data-processing problem.

## Fitting trees and forests

### Forests

**Wilcoxon random forest (WRF)**

``` r
# Wilcoxon Random Forest (WRF)
wrf <- myRF(y ~ ., data = train,
            ntree = 1000, 
            mtry = 2, # candidate predictors per split (default 2)
            nodesize = 10, # minimum leaf size (default 5)
            subsample_prop = 1 - exp(-1),# ~0.632 subsampling (default)
            replace = FALSE, # subsampling, not bootstrap (default)
            split_rule = "Wilcoxon")
```

**Quantile regression forest comparator**

Within this implementation, a QRF-style estimator uses SSE splitting and the same weighted-ECDF prediction rule:

``` r
# Quantile Regression Forest (QRF) = SSE splitting + weighted ECDF
qrf <- myRF(y ~ ., data = train, ntree = 1000, mtry = 2,
            nodesize = 10, split_rule = "SSE")
```

**Pinball-loss forest**

The numeric pinball split search uses an approximate grid of candidate cut positions for speed.

``` r
# Pinball-loss forest (targets specific quantile levels; approximate split search)
# (not recommended for this paper but optional)
pbf <- myRF(y ~ ., data = train, ntree = 1000, mtry = 2, nodesize = 10,
            split_rule = "pinball", pinball_taus = c(0.1, 0.5, 0.9))
```

*Key arguments:*

-   `ntree`: number of trees per forest, default 1000

-   `mtry`: candidate predictors per split, default 2

-   `nodesize`: minimum leaf size; default 5; should increase as training size increases

-   `max_depth`: maximum depth of tree, default 1000

-   `max_nodes`: maximum number of nodes, optional

-   `replace`: bootstrap if `TRUE`, subsample if `FALSE`; default `FALSE`

-   `subsample_prop`: subsampling proportion; default `1-exp(-1)`; required to be (0,1], recommended to be [0.5, 1]

-   `pinball_taus`: used only by the pinball-loss splitting rule, quantile(s) considered for pinball-loss impurity reduction.

### Single tree

A single tree is fit by setting `ntree = 1`, using the full sample `subsample_prop = 1` and `replace = FALSE`, and considering all predictors at each split `mtry = p`:

``` r
p <- ncol(train) - 1 # or the number of covariates

wilcoxon_tree <- myRF(y ~ ., data = train,
                      ntree = 1, # fit 1 tree only
                      subsample_prop = 1, replace = FALSE, # use the full sample
                      mtry = p, # consider all covariates
                      nodesize = 60,
                      split_rule = "Wilcoxon") # similar for other splitting rules
```
Although in the manuscript, we described a *pruning* method for the Wilcoxon regression tree, we have not yet implemented it in `R` as of September 25, 2026. We will implement it in the near future.


## Prediction

### Quantiles, means, medians, probabilities, CDFs, and weights

To predict conditional quantiles:

``` r
# Conditional quantiles
predict(wrf, # fitted WRF object
        newdata = test, # data frame for test observations
        what = "quantiles", # what you want to predict
        probs = c(0.1, 0.5, 0.9), # specific quantile levels
        linear_intpl = FALSE) # no linear interpolation
```
`linear_intpl = FALSE` by default and the predicted $\tau$th quantile follows: $\hat{Q}_\tau(x)= \inf \{y: \hat{F} (y \mid x) \geq \tau\}$. When `linear_intpl = TRUE`, the predicted quantile is obtained via linear interpolation, which is the same as `type = 4` for the `quantile()` function.

To predict conditional median, we can also use:

``` r
predict(wrf, newdata = test, what = "median")
```

To predict conditional mean:

``` r
# Conditional mean
predict(wrf, newdata = test, what = "mean")
```

To predict threshold probability:

``` r
# P(Y <= c | X) at several thresholds
predict(wrf, newdata = test, 
        what = "probability", # as we want to predict probabilities
        outcomes = c(200, 350, 500), # threshold value(s)
        lower.tail = TRUE, # TRUE for less than, FALSE for greater than
        include = TRUE) # TRUE for equal to, FALSE for strict relation
```

To predict $P(Y < c \mid X)$, set `lower.tail = TRUE` and `include = FALSE`.
To predict $P(Y \geq c \mid X)$, set `lower.tail = FALSE` and `include = TRUE`.
To predict $P(Y > c \mid X)$, set `lower.tail = FALSE` and `include = FALSE`.

To predict the conditional CDF or forest weights for a *single* observation: 

``` r
# Predicted CDF or forest weights for one profile
predict(wrf, newdata = test[1, ], what = "ECDF") # for ECDF
predict(wrf, newdata = test[1, ], what = "weights") # for forest weights
```



### Distributional evaluation

The following calls require the observed test outcomes:

``` r
# CRPS
crps <- predict(wrf, newdata = test, what = "CRPS",
                y_obsvd = y_test) # the observed outcome values for the test data

# Cramér-von Mises statistic for calibration
cvm <- predict(wrf, newdata = test, what = "CvM",
               y_obsvd = y_test)
```

`CRPS` returns one score per observation. `CvM` returns an aggregate Cramér–von Mises statistic. 

### Out-of-bag prediction

``` r
predict_oob.myRF(wrf, what = "quantiles",
                 probs = c(0.1, 0.5, 0.9), 
                 response_members = "all")

predict_oob.myRF(wrf, what = "probability",
                 outcomes = c(200, 350, 500),
                 response_members = "all")

predict_oob.myRF(wrf, what = "CRPS", response_members = "all")
```

`response_members` controls which observations are used to construct the terminal-node outcome distribution for out-of-bag prediction. With `response_members = "all"`, all training observations in the leaf contribute except the target observation itself; with `response_members = "inbag"`, only observations used to grow that tree contribute.

All OOB analyses reported in our manuscript use `response_members = "all"`.

## CRPS-based permutation variable importance

``` r
vi <- variable_importance_crps.myRF(
  wrf,
  R = 10, # number of repetitions
  seed = 2026, # randomization seed
  response_members = "all", # for OOB prediction
  sort = TRUE, # # if TRUE then output summary is ordered from largest to smallest importance
  verbose = TRUE) # if TRUE then prints progress messages

vi$summary
```

The fitted trees remain fixed. Within each tree, one predictor is permuted among out-of-bag observations, those observations are rerouted through the tree, and the increase in mean tree-level out-of-bag CRPS is calculated.

An out-of-bag pinball-loss diagnostic is also available:

``` r
oob_pinball_loss.myRF(wrf, tau = 0.5)
```

## CPM refinement

`myCPMForest()` fits an `rms::orm()` cumulative probability model (CPM) using terminal-node indicators from an already fitted `myRF` object. Please make sure the `rms` package has version >= 7.0.

A CPM-refined tree must be based on a fitted Wilcoxon regression tree structure. To

``` r
source("R/CPMtree.R")
cpm_tree <- myCPMForest(wilcoxon_tree, # the fitted Wilcoxon regression tree
                        family = "probit") # link function
```
`family` is the same link function as in the `orm.fit()` function in the `rms` package, and can be `logistic`, `probit`, `loglog`, `cloglog` or `cauchit`.

Trees that cannot support a CPM fit are skipped. Note that CPM refinement can be substantially slower than ECDF-based prediction.

For predictions:

``` r
predict(cpm_tree, newdata = test, what = "quantiles",
        probs = c(0.1, 0.5, 0.9))
predict(cpm_tree, newdata = test, what = "probability",
        outcomes = c(200, 350, 500))
predict(cpm_tree, newdata = test, what = "CRPS",
        y_obsvd = y_test)
```


## Toy example and sanity check

``` r
set.seed(1)
dat_train <- data.frame(
  y = rnorm(100),
  x1 = rnorm(100),
  x2 = factor(sample(c("a", "b"), 100, replace = TRUE)),
  x3 = rexp(100))
dat_test <- data.frame(y = 1, x1 = -1, x2 = "a", x3 = 2)

fit <- myRF(y ~ ., dat_train, ntree = 20, mtry = 2,
            nodesize = 5, split_rule = "Wilcoxon")

predict(fit, dat_test, what = "quantiles")
#       q_0.1      q_0.5    q_0.9
# 1 -1.044135 0.07456498 1.160403

pred_ECDF <- predict(fit, dat_test[1, ], what = "ECDF")
head(pred_ECDF)
#     y.sorted       ECDF
# 1  -2.214700 0.01041667
# 4  -1.523567 0.03162879
# 6  -1.377060 0.05859307
# 8  -1.253633 0.07039863
# 10 -1.129363 0.08379149
# 11 -1.044135 0.10611291

```

## Reference

Shi, D., Shepherd, B.E., Li, C., (2026) Wilcoxon Random Forests for Robust Distributional Prediction.

Harrell Jr FE (2025). _rms: Regression Modeling Strategies_. doi:10.32614/CRAN.package.rms <https://doi.org/10.32614/CRAN.package.rms>, R package version 8.1-0, <https://CRAN.R-project.org/package=rms>.

Bates D, Maechler M, Jagan M (2025). _Matrix: Sparse and Dense Matrix Classes and Methods_. doi:10.32614/CRAN.package.Matrix <https://doi.org/10.32614/CRAN.package.Matrix>, R package version 1.7-4, <https://CRAN.R-project.org/package=Matrix>.
