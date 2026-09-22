# Cumulative Probability Model (CPM) refinement for fitted WRF/QRF-style trees
# This file fits a cumulative probability model inside each tree of a fitted
# myRF object.  The CPM uses terminal-node indicators as predictors and then
# averages the resulting tree-level CDFs across trees.
#
# Source the myRF code before this file, because the CPM code uses the same
# vectorized tree-traversal helpers used by myRF prediction.

.check_rms_ver <- function() {
  if (!requireNamespace("rms", quietly = TRUE)) {
    stop("Package 'rms' is required (>= 7.0). Install it with install.packages('rms').")
  }
  if (utils::packageVersion("rms") < "7.0-0") {
    stop("Package 'rms' must be version 7.0 or newer.")
  }
}

# Build dummy variables for terminal-node IDs.
#
# The first terminal node is used as the reference level.  At prediction time,
# reference_levels ensures that the new dummy matrix has the same columns as the
# matrix used when fitting the tree-level CPM.
.make_cpm_dummies <- function(leaf_id_vec, reference_levels = NULL) {
  if (is.null(reference_levels)) {
    lf <- factor(leaf_id_vec)
  } else {
    lf <- factor(leaf_id_vec, levels = reference_levels)
  }
  
  if (length(levels(lf)) < 2) {
    return(list(Xd = data.frame(row.names = seq_along(leaf_id_vec)),
                levels = levels(lf)))
  }
  
  MM <- model.matrix(~ lf)
  MM <- MM[, -1, drop = FALSE] # Drop intercept to avoid collinearity.
  colnames(MM) <- paste0("Node_", levels(lf)[-1])
  
  list(Xd = as.data.frame(MM), levels = levels(lf))
}

# Fit a CPM forest from an already fitted myRF object.
#
# Each base tree is kept fixed. The tree partition was selected using that tree's
# subsample/bootstrap sample. Following the current forest estimator, all original
# training observations are then routed through the fitted tree, terminal-node
# indicators are created, and an rms::orm() cumulative probability model is fit.
# Trees that cannot support a CPM fit are skipped.
myCPMForest <- function(base_forest,
                        family = c("logistic", "probit", "cloglog", "loglog", "cauchit"),
                        verbose = TRUE) {
  if (!inherits(base_forest, "myRF")) stop("base_forest must be a myRF object.")
  if (!exists(".get_leaf_ids_vectorized", mode = "function")) {
    stop("Source the myRF code before fitting a CPM forest.")
  }

  .check_rms_ver()
  family <- match.arg(family)

  y <- base_forest$y
  X <- base_forest$X
  ntree <- base_forest$ntree
  RF <- base_forest$RF

  if (length(unique(y)) < 2L) {
    stop("CPM refinement requires at least two unique outcome values.")
  }

  cpm_forest <- vector("list", ntree)
  kept <- 0L

  for (t in seq_len(ntree)) {
    tree <- RF[[t]]

    # Route all original training observations through the fixed tree. The tree
    # structure itself was constructed only from that tree's in-bag observations.
    terminal_nodes_all <- if (!is.null(tree$leaf_id_by_obs)) {
      tree$leaf_id_by_obs
    } else {
      .get_leaf_ids_vectorized(tree, X)
    }

    if (length(unique(terminal_nodes_all)) < 2L) {
      if (verbose) message(sprintf("[Tree %d] Skipped: only one terminal node among training observations.", t))
      next
    }

    # Encode terminal-node membership as dummy variables.
    dd <- .make_cpm_dummies(terminal_nodes_all)
    dummy_vars <- colnames(dd$Xd)
    if (length(dummy_vars) == 0L) next

    cpm_data <- data.frame(Y = y, dd$Xd, check.names = FALSE)
    cpm_formula <- as.formula(paste("Y ~", paste(dummy_vars, collapse = " + ")))

    # Fit the tree-level cumulative probability model.
    fit <- try(
      rms::orm(cpm_formula, data = cpm_data, family = family, x = TRUE, y = TRUE),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || isTRUE(fit$fail)) {
      if (verbose) message(sprintf("[Tree %d] orm() failed; skipped.", t))
      next
    }

    kept <- kept + 1L
    cpm_forest[[kept]] <- list(
      tree = tree,
      cpm_fit = fit,
      terminal_node_levels = dd$levels,
      tree_id = t,
      y_support = sort(unique(y))
    )
  }

  if (kept == 0L) stop("No trees were successfully fitted with CPM.")
  if (verbose) message(sprintf("CPM fits kept: %d/%d trees.", kept, ntree))

  cpm_forest <- cpm_forest[seq_len(kept)]

  structure(list(
    cpm_forest = cpm_forest,
    base_forest = base_forest,
    ntree = kept,
    family = family,
    y = y,
    X = X,
    yname = base_forest$yname,
    xnames = base_forest$xnames,
    unique_y = sort(unique(y))
  ), class = "myCPMForest")
}

# CRPS for predictions represented by a step CDF on y_grid.
.calc_crps_cpm <- function(F_mat, y_grid, y_obsvd) {
  if (length(y_obsvd) != nrow(F_mat)) {
    stop("Length of 'y_obsvd' must match the number of rows in the prediction matrix.")
  }
  
  K <- length(y_grid)
  
  # Convert the CDF matrix to probability masses on the support points.
  if (K > 1L) {
    W_mat <- cbind(F_mat[, 1L],
                   F_mat[, 2:K, drop = FALSE] - F_mat[, 1:(K - 1L), drop = FALSE])
  } else {
    W_mat <- F_mat
  }
  
  # Equivalent CRPS representation: E|X - y_obs| - 0.5 E|X - X'|.
  dist_mat <- abs(outer(y_obsvd, y_grid, "-"))
  term1 <- rowSums(W_mat * dist_mat)
  
  if (K > 1L) {
    dy <- diff(y_grid)
    term2 <- as.vector((F_mat[, -K, drop = FALSE] *
                          (1 - F_mat[, -K, drop = FALSE])) %*% dy)
  } else {
    term2 <- 0
  }
  
  term1 - term2
}

# Cramer-von Mises statistic from PIT values, or probability-scale residual (PSR) values, for a step CDF.
.calc_cvm_cpm <- function(F_mat, y_grid, y_obsvd, PSR) {
  if (length(y_obsvd) != nrow(F_mat)) {
    stop("Length of 'y_obsvd' must match the number of rows in the prediction matrix.")
  }
  
  n_test <- length(y_obsvd)
  idx_le <- findInterval(y_obsvd, y_grid) # P(Y <= y_obs)
  idx_lt <- findInterval(y_obsvd, y_grid, left.open = TRUE)        # P(Y < y_obs)
  
  p_le <- numeric(n_test)
  valid_le <- idx_le > 0L
  if (any(valid_le)) p_le[valid_le] <- F_mat[cbind(which(valid_le), idx_le[valid_le])]
  
  p_lt <- numeric(n_test)
  valid_lt <- idx_lt > 0L
  if (any(valid_lt)) p_lt[valid_lt] <- F_mat[cbind(which(valid_lt), idx_lt[valid_lt])]
  
  # Mid-PIT: F(y-) + 0.5 * P(Y = y).  This is useful for discrete or tied outcomes.
  U_vec <- p_lt + 0.5 * (p_le - p_lt)
  
  if (PSR) {
    return(2 * U_vec - 1)
  }
  
  U_sorted <- sort(U_vec)
  i <- seq_len(n_test)
  sum((U_sorted - (2 * i - 1) / (2 * n_test))^2) + 1 / (12 * n_test)
}

# Evaluate lower- or upper-tail probabilities from a step-CDF matrix.
.cpm_prob_from_ecdf_matrix <- function(F_mat, y_grid, outcomes,
                                       lower.tail = TRUE,
                                       include = TRUE,
                                       linear_intpl = FALSE) {
  n_new <- nrow(F_mat)
  P <- matrix(NA_real_, nrow = n_new, ncol = length(outcomes))
  
  if (linear_intpl) {
    if (length(y_grid) == 1L) {
      p_i <- as.numeric(outcomes >= y_grid[1L])
      P[,] <- rep(if (lower.tail) p_i else 1 - p_i, each = n_new)
      colnames(P) <- as.character(outcomes)
      return(P)
    }

    for (i in seq_len(n_new)) {
      p_i <- approx(x = y_grid, y = F_mat[i, ], xout = outcomes,
                    method = "linear", rule = 2, ties = "ordered")$y
      p_i[outcomes < y_grid[1L]] <- 0
      P[i, ] <- if (lower.tail) p_i else 1 - p_i
    }
    colnames(P) <- as.character(outcomes)
    return(P)
  }
  
  idx_le <- findInterval(outcomes, y_grid) # # {y <= c}
  idx_lt <- findInterval(outcomes, y_grid, left.open = TRUE)        # # {y <  c}
  
  for (k in seq_along(outcomes)) {
    F_le <- if (idx_le[k] > 0L) F_mat[, idx_le[k]] else rep(0, n_new)
    F_lt <- if (idx_lt[k] > 0L) F_mat[, idx_lt[k]] else rep(0, n_new)
    
    if (lower.tail) {
      P[, k] <- if (include) F_le else F_lt
    } else {
      P[, k] <- if (include) 1 - F_lt else 1 - F_le
    }
  }
  
  colnames(P) <- as.character(outcomes)
  P
}

# Predict from a fitted CPM forest.
#
# The CPM forest returns the same main prediction types as myRF: ECDFs,
# quantiles, means, threshold probabilities, CRPS, PSR, and CvM.
predict.myCPMForest <- function(object, newdata,
                                what = c("ECDF", "quantiles", "mean", "probability", "CRPS", "PSR", "CvM"),
                                probs = c(0.1, 0.5, 0.9),
                                outcomes = NULL,
                                y_obsvd = NULL,
                                lower.tail = TRUE,
                                include = TRUE,
                                linear_intpl = FALSE,
                                ...) {
  what <- match.arg(what)
  newdata <- as.data.frame(newdata)[, object$xnames, drop = FALSE]
  if (anyNA(newdata)) stop("newdata contains NA(s).")
  
  n_new <- nrow(newdata)
  unique_y <- object$unique_y
  n_unique <- length(unique_y)
  
  ecdf_matrix <- matrix(0, nrow = n_new, ncol = n_unique)
  tree_contributions <- 0L
  
  for (tree_obj in object$cpm_forest) {
    # Route new observations through the fixed base tree.
    terminal_nodes <- .get_leaf_ids_vectorized(tree_obj$tree, newdata)
    
    # Build new-data dummies with the same terminal-node levels used for fitting.
    new_dummies <- .make_cpm_dummies(terminal_nodes, tree_obj$terminal_node_levels)$Xd
    if (ncol(new_dummies) == 0L && length(tree_obj$terminal_node_levels) >= 2L) {
      new_dummies <- as.data.frame(
        matrix(0, nrow = n_new, ncol = length(tree_obj$terminal_node_levels) - 1L)
      )
      colnames(new_dummies) <- paste0("Node_", tree_obj$terminal_node_levels[-1])
    }
    
    # Predict the tree-level CPM.  For rms::orm(), type = "fitted" returns
    # exceedance probabilities on the model's outcome support; the conversion
    # below gives P(Y <= y) at each support point.
    p_ge <- try(predict(tree_obj$cpm_fit, newdata = new_dummies, type = "fitted"),
                silent = TRUE)
    if (inherits(p_ge, "try-error")) next
    
    yt <- tree_obj$cpm_fit$yunique
    n_tree_unique <- length(yt)
    
    if (is.vector(p_ge)) {
      if (n_new == 1L) {
        p_ge <- matrix(p_ge, nrow = 1L)
      } else {
        p_ge <- matrix(p_ge, ncol = 1L)
      }
    } else {
      p_ge <- as.matrix(p_ge)
    }
    
    p_le_tree <- matrix(1, nrow = n_new, ncol = n_tree_unique)
    if (n_tree_unique >= 2L) {
      p_le_tree[, 1:(n_tree_unique - 1L)] <- 1 - p_ge
    }
    
    # Map this tree's CDF to the global outcome grid.  For global grid values
    # below this tree's smallest support point, the tree-level CDF is zero.
    idx_map <- findInterval(unique_y, yt)
    valid_map <- idx_map > 0L
    if (any(valid_map)) {
      ecdf_matrix[, valid_map] <-
        ecdf_matrix[, valid_map, drop = FALSE] +
        p_le_tree[, idx_map[valid_map], drop = FALSE]
    }

    tree_contributions <- tree_contributions + 1L
  }
  
  if (tree_contributions == 0L) stop("No trees produced valid predictions.")
  
  # Average tree-level CDFs and enforce valid monotone CDF values.
  final_ecdf <- ecdf_matrix / tree_contributions
  final_ecdf <- t(apply(final_ecdf, 1, function(z) {
    z <- pmin(pmax(z, 0), 1)
    if (length(z) >= 2L) {
      for (j in 2:length(z)) if (z[j] < z[j - 1L]) z[j] <- z[j - 1L]
    }
    z
  }))
  
  if (what == "ECDF") {
    if (n_new != 1L) stop("ECDF output is only supported for a single prediction row.")
    return(data.frame(y = unique_y, ECDF = final_ecdf[1, ]))
  }
  
  if (what == "quantiles") {
    Q0 <- apply(final_ecdf, 1, function(Fi) {
      sapply(probs, function(p) {
        if (p <= Fi[1L]) return(unique_y[1L])
        if (p >= Fi[n_unique]) return(unique_y[n_unique])

        idx <- which(Fi >= p)[1L]
        if (is.na(idx)) return(unique_y[n_unique])

        if (linear_intpl && idx > 1L) {
          y1 <- unique_y[idx - 1L]
          y2 <- unique_y[idx]
          f1 <- Fi[idx - 1L]
          f2 <- Fi[idx]
          if (f2 > f1) {
            return(y1 + (y2 - y1) * (p - f1) / (f2 - f1))
          }
        }
        unique_y[idx]
      })
    })

    if (length(probs) == 1L) {
      Q <- matrix(Q0, ncol = 1L)
    } else {
      Q <- t(Q0)
    }
    colnames(Q) <- paste0("q_", probs)
    return(as.data.frame(Q))
  }
  if (what == "mean") {
    means <- apply(final_ecdf, 1, function(Fi) {
      w <- c(Fi[1L], diff(Fi))
      w[w < 0] <- 0
      sw <- sum(w)
      if (sw > 0) w <- w / sw
      sum(unique_y * w)
    })
    return(means)
  }
  
  if (what == "probability") {
    if (is.null(outcomes)) stop("Provide 'outcomes' when what = 'probability'.")
    P <- .cpm_prob_from_ecdf_matrix(final_ecdf, unique_y, outcomes,
                                    lower.tail = lower.tail,
                                    include = include,
                                    linear_intpl = linear_intpl)
    return(as.data.frame(P))
  }
  
  if (what == "CRPS") {
    if (is.null(y_obsvd)) stop("You must provide 'y_obsvd' (observed y) for CRPS.")
    return(.calc_crps_cpm(final_ecdf, unique_y, y_obsvd))
  }
  
  if (what == "CvM" || what == "PSR") {
    if (is.null(y_obsvd)) stop("You must provide 'y_obsvd' (observed y) for CvM or PSR.")
    return(.calc_cvm_cpm(final_ecdf, unique_y, y_obsvd, PSR = (what == "PSR")))
  }
}
