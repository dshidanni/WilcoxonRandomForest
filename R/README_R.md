# Wilcoxon Random Forest (WRF) and CPM-Refined Regression Trees

*Author: me!*

These R codes fit regression trees and forests under three splitting rules: 1. sum of squared errors (SSE), 2. rank-based Wilcoxon, 3. pinball-loss. The R codes also estimate conditional distributions via forest-weighed empirical CDFs, and provides a cumulative probability model (CPM) refined of the fitted Wilcoxon regression trees. 

---

## 1. Contents

| File        | What it provides                                                                                                                                                                                                                                                                                                        |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `myRF.R`    | Core code: fit trees/forests (`myRF`), predict (`predict.myRF`), out-of-bag prediction (`predict_oob.myRF`), CRPS-based permutation variable importance (`variable_importance_crps.myRF`), OOB pinball-loss summary (`oob_pinball_loss.myRF`), and the distributional evaluation metrics (CRPS, Cramér–von Mises, PSR). |
| `CPMtree.R` | CPM refinement of a fitted forest (`myCPMForest`) and its prediction method (`predict.myCPMForest`). Source **after** `myRF.R`.                                                                                                                                                                                         |

> Notes: `predict.myRF` and `predict.myCPMForest` are S3 methods, so you can
> call `predict(object, ...)`. Other functions (`predict_oob.myRF`,
> `variable_importance_crps.myRF`, `oob_pinball_loss.myRF`) are **not** S3 generics, so please call them by their full name.

---

## 2. Requirements

- **R** (developed under R 4.5.3; the paper used the version cited as `r2026`).
- **Matrix** (sparse forest-weight matrices). The paper used version **1.7-4**.
- **rms** (only needed for the cumulative probability model (CPM) to refine Wilcoxon regression tree). This paper used version 8.1-1, and requires version $\geq$ **7.0-0**.

```r
install.packages(c("Matrix", "rms"))
```

```r
> sessionInfo()
R version 4.5.3 (2026-03-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] rms_8.1-0   Hmisc_5.2-5

loaded via a namespace (and not attached):
 [1] sandwich_3.1-1      generics_0.1.4      stringi_1.8.7       lattice_0.22-9      digest_0.6.39      
 [6] magrittr_2.0.4      evaluate_1.0.5      grid_4.5.3          RColorBrewer_1.1-3  mvtnorm_1.3-3      
[11] fastmap_1.2.0       Matrix_1.7-4        nnet_7.3-20         backports_1.5.0     Formula_1.2-5      
[16] survival_3.8-6      multcomp_1.4-29     gridExtra_2.3       scales_1.4.0        TH.data_1.1-5      
[21] codetools_0.2-20    cli_3.6.5           rlang_1.1.7         splines_4.5.3       base64enc_0.1-6    
[26] otel_0.2.0          tools_4.5.3         MatrixModels_0.5-4  SparseM_1.84-2      checkmate_2.3.4    
[31] htmlTable_2.4.3     dplyr_1.2.0         colorspace_2.1-2    ggplot2_4.0.2       vctrs_0.7.1        
[36] R6_2.6.1            rpart_4.1.24        zoo_1.8-15          polspline_1.1.25    lifecycle_1.0.5    
[41] stringr_1.6.0       htmlwidgets_1.6.4   MASS_7.3-65         foreign_0.8-91      cluster_2.1.8.2    
[46] pkgconfig_2.0.3     pillar_1.11.1       gtable_0.3.6        glue_1.8.0          data.table_1.18.2.1
[51] xfun_0.56           tibble_3.3.1        tidyselect_1.2.1    rstudioapi_0.18.0   knitr_1.51         
[56] farver_2.1.2        nlme_3.1-168        htmltools_0.5.9     rmarkdown_2.30      compiler_4.5.3     
[61] quantreg_6.1        S7_0.2.1
```

---

## 3. Loading the code

```r
source("myRF.R") # always
source("CPMtree.R")  # only if you use the CPM refinement (requires rms>=7.0)
```

---

## 4. Data requirements

- **Response variable must be numeric**: continuous, count, or an ordered/discrete outcome
  that has been **coded to numeric values preserving its order**. Convert ordered
  factors to numeric *before* calling `myRF()`; the function stops on a non-numeric
  response. (Detection-limit / mixed outcomes should be coded so the order is
  respected, e.g. a below-limit code strictly below all measured values.)
- **Covariates** can be numeric or categorical (factor / character). Numeric splits
  use `x <= s`; categorical splits store the set of levels sent left. A category not
  seen during training is routed to the right child during prediction.
- **Missing values are not imputed.** `myRF()` stops if the model frame contains `NA`.
  `Inf`/`-Inf` stop SSE and pinball fits and warn (but proceed) for Wilcoxon.

---

## 5. Leaf-population convention (please read)

For prediction on **new data**, each terminal node's outcome distribution is formed
from **all** training observations that fall into that node (every training point is
dropped down every fitted tree). This is the **all-sample** estimator and matches the
QRF convention of Meinshausen (2006): trees are grown on subsamples, but the leaf ECDFs
use all *n* observations.

**Out-of-bag** prediction (`predict_oob.myRF`, and the permutation variable importance
built on it) instead forms each node's distribution from that tree's **in-bag**
subsample only (`response_members = "inbag"`, the default). This is deliberate and is a
*different* estimator: for an OOB observation, using all training observations would
leak that observation's own outcome into its own predicted distribution.

---

## 6. Fitting

### 6.1 Forests: WRF, QRF, and pinball

```r
# Wilcoxon Random Forest (WRF)
wrf <- myRF(y ~ ., data = train,
            ntree = 1000, 
            mtry = 2, # candidate predictors per split (default 2)
            nodesize = 10, # minimum leaf size (default 5)
            subsample_prop = 1 - exp(-1),# ~0.632 subsampling (default)
            replace = FALSE, # subsampling, not bootstrap (default)
            split_rule = "Wilcoxon")

# Quantile Regression Forest (QRF) = SSE splitting + weighted ECDF
qrf <- myRF(y ~ ., data = train, ntree = 1000, mtry = 2,
            nodesize = 10, split_rule = "SSE")

# Pinball-loss forest (targets specific quantile levels; approximate split search)
# (not recommended for this paper but optional)
pbf <- myRF(y ~ ., data = train, ntree = 1000, mtry = 2, nodesize = 10,
            split_rule = "pinball", pinball_taus = c(0.1, 0.5, 0.9))
```

Key arguments: 
* `ntree`: number of trees per forest, default 1000
* `mtry`: candidate predictors per split, default 2
* `nodesize`: minimum leaf size; default 5; should increase as training size increases
* `max_depth`: maximum depth of tree, default 100
* `max_nodes`: maximum number of nodes, optional
* `replace`: bootstrap if `TRUE`, subsample if `FALSE`; default `FALSE`
* `subsample_prop`: subsampling proportion; default `1-exp(-1)`; required to be (0,1], recommended to be $[0.5, 1]$
* `pinball_taus`: quantile(s) considered for pinball-loss impurity reduction, used only by the pinball-loss splitting rule.


### 6.2 Single trees

Single trees are forests with one tree grown on the full training sample and all
predictors considered at each split:

```r
tree_wilcoxon <- myRF(y ~ ., data = train, ntree = 1,
                      mtry = <p>, # all predictors (p = number of covariates)
                      nodesize = 60,
                      subsample_prop = 1,  # use the full sample (no subsampling)
                      replace = FALSE,
                      split_rule = "Wilcoxon")
tree_cart  <- myRF(..., split_rule = "SSE") # CART regression tree
tree_pinbl <- myRF(..., split_rule = "pinball")  # pinball-loss tree
```


### 6.3 CPM-refined regression tree

Fit a cumulative probability model inside each fitted Wilcoxon regression tree, using terminal-node
indicators as predictors, then average the tree-level CDFs:

```r
cpm_probit <- myCPMForest(tree_wilcoxon, 
						  family = "probit") # link function
```

`myCPMForest` takes a fitted Wilcoxon regression tree, (a `myRF` object) and does not change the tree structure. `myCPMForest` can also take a fitted WRF, although it is not recommended for this paper. Trees that cannot support a CPM fit (e.g. only one terminal node among in-bag observations) are skipped.
`family` argument is the same as the `family` argument in `rms::orm()` function: it specifies the family of the link function and takes character values of `"logistic", "probit", "loglog", "cloglog", "cauchit"`. See `rms` package manual for details.

---

## 7. Prediction

### 7.1 From a forest (`predict.myRF`)

```r
predict(wrf, newdata, what = "quantiles", probs = c(0.1, 0.5, 0.9))
predict(wrf, newdata, what = "mean")
predict(wrf, newdata, what = "median")
predict(wrf, newdata, what = "probability", outcomes = c(0, 200, 350))  # P(Y <= c) by default
predict(wrf, newdata, what = "ECDF")     # nrow(newdata) == 1 only
predict(wrf, newdata, what = "weights")  # nrow(newdata) == 1 only
```

Options:
- `probs` — quantile levels for `what = "quantiles"`/`"median"`.
- `outcomes` — thresholds for `what = "probability"`.
- `lower.tail`, `include` — control which tail and whether the boundary is included:
  `lower.tail = TRUE, include = TRUE` → P(Y ≤ c); `include = FALSE` → P(Y < c);
  `lower.tail = FALSE` → upper-tail counterparts.
- `linear_intpl` — `FALSE` (default) gives step-function quantiles/probabilities;
  `TRUE` linearly interpolates between support points.

Distributional scoring (require the observed outcomes via `y_obsvd`):

```r
predict(wrf, newdata, what = "CRPS", y_obsvd = y_test)  # continuous ranked probability score
predict(wrf, newdata, what = "CvM",  y_obsvd = y_test)  # Cramér–von Mises PIT statistic
predict(wrf, newdata, what = "PSR",  y_obsvd = y_test)  # probability-scale residuals
predict(wrf, newdata, what = "PSR discrete", y_obsvd = y_test)
```

### 7.2 From a CPM-refined tree/forest (`predict.myCPMForest`)

```r
predict(cpm_probit, newdata, what = "ECDF") # newdata must have single row
predict(cpm_probit, newdata, what = "quantiles", probs = c(.1,.5,.9))
predict(cpm_probit, newdata, what = "mean")
predict(cpm_probit, newdata, what = "probability", outcomes = c(200,350,500))
predict(cpm_probit, newdata, what = "CRPS", y_obsvd = y_test)
predict(cpm_probit, newdata, what = c("CvM","PSR")[1], y_obsvd = y_test)
```

---

## 8. Out-of-bag prediction

For each subject, only trees in which that subject was out-of-bag contribute; the
terminal-node distribution uses in-bag members by default (`response_members = "inbag"`):

```r
predict_oob.myRF(wrf, what = "quantiles", probs = c(.1,.5,.9))
predict_oob.myRF(wrf, what = "mean")
predict_oob.myRF(wrf, what = "probability", outcomes = c(0,1))
predict_oob.myRF(wrf, what = "CRPS")   # uses the training outcomes as y_obsvd internally
```

---

## 9. Variable importance

CRPS-based OOB permutation importance (percentage increase in mean OOB CRPS when a
covariate is permuted among OOB observations; trees are kept fixed):

```r
vi <- variable_importance_crps.myRF(wrf,
                                    R = 10,        # permutation repetitions per covariate
                                    seed = 1,      # reproducibility
                                    sort = TRUE)
vi$summary       # variable, baseline/permuted mean CRPS, VI_CRPS, VI_percent
```

An OOB mean pinball-loss diagnostic is also available:

```r
oob_pinball_loss.myRF(wrf, tau = 0.5)
```

---

## 10. Reproducing the paper's results

Settings used (fill in any exact values you changed):

- **Forests:** `ntree = 1000`, `mtry = 2`, `subsample_prop = 1 - exp(-1)`,
  `replace = FALSE`. `nodesize`: Scenario 1 uses `{10, 20, 50}` and Scenario 2 uses
  `{5, 10, 25}`, matching `n_train ∈ {1000, 2000, 5000}`.
- **Single trees:** `ntree = 1`, all predictors per split, `nodesize ∈ {60, 120, 300}`.
- **QRF vs WRF:** `split_rule = "SSE"` vs `"Wilcoxon"`.
- **CPM trees/forests:** `family = "probit"` (correct) and `"logistic"` (misspecified).
- **Evaluation:** a fixed test set of `n_test = 1000`; point/interval and distributional
  metrics via `predict(..., what = ...)`; forest comparisons averaged over 250 replicates,
  single-tree comparisons over 1000.
- **Reproducibility:** set a seed (`set.seed(...)`) at the start of each replicate, and
  pass `seed = ` to `variable_importance_crps.myRF`.

### Not included in these two files (you must add them for a full replication)

- **Simulation driver scripts** — the data-generating processes (Scenarios 1.1, 1.2, 2),
  the replicate loop, and the code that assembles the RMSE/CRPS/CvM tables are **not**
  in `myRF.R` / `CPMtree.R`. Add your driver script(s) to the repository.
- **Application data and scripts** — the clinical data are not distributed here; include a
  data-availability statement and the analysis script.
- **Cost-complexity pruning** — the rank-based pruning described in the supplement is
  **not** implemented in these files (tree size is controlled by `nodesize` / `max_depth`
  / `max_nodes`). If any reported result used pruning, include that code as well.

---

## 11. Notes on this release

Relative to the originally submitted scripts, this version makes **no change to any
reported numerical result**. The edits are: (i) require a numeric response instead of
silently coercing it; (ii) fix the linear-interpolation quantile helper so it ignores
zero-weight support points (this path is only used when `linear_intpl = TRUE`, which the
reported analyses did not use); (iii) performance-only changes that leave outputs
bit-identical — iterate sparse weight matrices by row instead of densifying them, and use
row-compressed storage for the per-row evaluation loops; (iv) a vectorized (identical)
CDF monotonization in the CPM predictor; and (v) added explanatory comments and a guard
on the `rms::orm` output layout. The step-function predictions, CRPS, Cramér–von Mises,
PSR, means, OOB predictions, and variable importance are unchanged.
