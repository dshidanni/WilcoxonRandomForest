# Simulation studies

This directory includes the R scripts for the simulation studies in Section 4 and Supplementary Section S2 of **“Wilcoxon Random Forests for Robust Distributional Prediction.”** The studies compare Wilcoxon regression trees and forests with SSE-based and pinball-loss-based alternatives for conditional quantile, threshold-probability, mean, and full-distribution prediction.

All scripts assume that the working directory is the root directory of this repository. File paths are specified relative to the repository root.

## Existing scripts

| File | Purpose |
|------------------------------------|------------------------------------|
| `simulation_scenario1.R` | Generates and analyzes Scenarios 1.1 and 1.2, including single-tree comparisons, QRF-versus-WRF forest comparisons, and CPM-refined WRF. |
| `simulation_scenario2.R` | Generates and analyzes Scenario 2, a smooth nonlinear model with right-skewed heteroscedastic errors. |
| `simulation_results.R` | Reads simulation outputs, calculates performance summaries, and writes the tables used in the manuscript and supplement. |

## Simulation design

### Common predictors

Each scenario uses six mutually independent predictors:

-   $X_1 \sim N(0, 1)$;
-   $X_2 \sim \text{Bernoulli}(p=0.5)$;
-   $X_3 \sim t_{\mathrm{df} = 10}$;
-   $X_4 \sim \text{Uniform} (0, 1)$;
-   $X_5 \sim \text{Beta} (3, 2)$;
-   $X_6 \sim \text{Poisson} (1)$.

$X_1$ to $X_4$ affect the outcome. $X_5$ and $X_6$ are noise predictors.

### Scenario 1.1: tree-structured, symmetric outcome

The conditional mean is piecewise constant on seven terminal regions, each with probability 1/7. The latent outcome is $$Y^* = \mu_1(\mathbf{X}) + \epsilon, \epsilon \sim N(0, 1).$$

Scenario 1.1 uses $Y = Y^*$. This scenario is symmetric and homoscedastic and was designed to be favorable to SSE-based splitting.

### Scenario 1.2: monotone transformation of Scenario 1.1

Scenario 1.2 uses the strictly increasing transformation $$h(y) = F^{-1}_{\chi^2_{\mathrm{df} = 15}}{\Phi(y)}.$$

The transformed outcome is right-skewed and heteroscedastic, while preserving the outcome ranks. Under identical randomization and tuning parameters, the Wilcoxon tree/forest partitions should therefore be unchanged between Scenarios 1.1 and 1.2.

### Scenario 2: smooth nonlinear, skewed, heteroscedastic outcome

The outcome is $Y = \mu_2(\mathbf{X}) + \sigma( \mathbf{X}) \cdot \epsilon$, $\epsilon = Z - 3$ where $Z \sim \chi^2_{\text{df} = 3}$, with

$$\mu_2(\mathbf{X}) = 0.5 + X_1 + 0.8X_2 + 0.3X_3 + 0.7\sin(2 \pi X_4) - 0.5X_1 X_2, \sigma(X) = 0.4 + 0.3I(X_4 > 0.5) + 0.5I(X_1 X_3 > 0).$$

This setting assesses performance away from a truly tree-structured regression function.

## Methods and settings

### Training and test samples

-   Training sizes: `n_train = 1000, 2000, 5000`.
-   Fixed test covariate set: `n_test = 1000`.
-   Test outcomes are generated independently within each Monte Carlo replicate, conditional on the fixed test covariates.
-   Single-tree comparisons: 1000 replicates.
-   Forest comparisons: 250 replicates.

### Single-tree comparisons

The Scenario 1 single-tree study compares:

1.  SSE regression tree;
2.  Wilcoxon regression tree;
3.  pinball-loss tree targeting `tau = 0.1, 0.5, 0.9` (minimizing the sum of pinball losses);
4.  CPM-refined Wilcoxon tree with probit link (correctly specified);
5.  CPM-refined Wilcoxon tree with logistic link (misspecified).

Each tree uses the full training sample, all six predictors as split candidates, and minimum terminal-node sizes `60, 120, 300` for training sizes `1000, 2000, 5000`. All five methods use the same randomization seed.

The pinball implementation uses an approximate grid of numeric split positions.

### Forest comparisons

The primary forest comparison is:

-   QRF-style estimator: SSE splitting plus forest-weighted ECDF prediction;
-   WRF: Wilcoxon splitting plus forest-weighted ECDF prediction.

Each forest uses:

-   1000 trees;
-   subsampling without replacement with subsampling fraction `1 - exp(-1)`;
-   `mtry = 2`;
-   the same randomization seeds.

Minimum node sizes are:

| Scenario   | `n_train = 1000` | `n_train = 2000` | `n_train = 5000` |
|------------|-----------------:|-----------------:|-----------------:|
| Scenario 1 |               10 |               20 |               50 |
| Scenario 2 |                5 |               10 |               25 |

### CPM-refined WRF

This exploratory comparison used the first 100 Monte Carlo replicates from the Scenario 1 forest experiment, with `n_train = 2000`, `n_test = 1000`, 1000 trees, `mtry = 2`, a subsampling fraction of 0.632, and `n_node = 20`.

Within each replicate, the WRF and the CPM-refined WRF used identical Wilcoxon tree partitions and forest randomization. The ECDF-based WRF used the terminal-node ECDF for each tree, whereas the CPM-refined WRF fitted a probit-link CPM with terminal-node indicators and used its fitted conditional CDF as the tree-level prediction.

## Estimands and performance measures

For each test observation, the scripts estimate:

-   conditional quantiles at `tau = 0.1, 0.5, 0.9`;
-   conditional mean;
-   in Scenario 1, a threshold probability corresponding to latent threshold zero;
-   the full predictive distribution.

Summaries include:

-   RMSE and bias of conditional quantiles, means, and threshold probabilities;
-   pinball loss at each reported quantile;
-   Brier score for the Scenario 1 threshold probability;
-   mean and median CRPS;
-   Cramér–von Mises calibration statistic;
-   80% prediction-interval coverage and median width;
-   empirical non-exceedance rates at 0.1 and 0.9;
-   elapsed fitting time.

The simulation results are displayed in Tables 1-2 and Tables S1-S18.


## Computation details

All code was implemented in R 4.5.3 with use of the `rms` package for CPM-based approaches (version 8.1-1), and the `Matrix` package for matrix computations (version 1.7-4). Timing was performed on an Intel Core Ultra 7-265 processor with 32 GB RAM without parallelization.

