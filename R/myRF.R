# Wilcoxon Random Forest (WRF) and related (Quantile) Regression Forest utilities
# 
# This file contains the core code used to fit distributional random forests
# with three splitting rules: squared-error (SSE), Wilcoxon rank-based, and
# pinball-loss splitting. 
# Prediction is based on forest weights and weighted
# empirical CDFs, following the quantile regression forest framework.
#
# Public user-facing functions:
#   myRF(): fit a forest
#   predict.myRF(): predict quantiles, probabilities, ECDFs, etc.
#   predict_oob.myRF(): out-of-bag predictions and diagnostics
#   variable_importance_crps.myRF(): CRPS-based OOB permutation importance
#   oob_pinball_loss.myRF(): OOB pinball-loss summary
#

# requires "Matrix" package
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required. Install it with install.packages('Matrix').")
}


# Fit a distributional random forest using one of three splitting rules.
myRF <- function(formula, # outcome ~ covariates
                 data, # training set
                 ntree=1000, # number of trees in forest
                 mtry=2, # number of candidate predictors sampled at each node,
                 # a common choice is sqrt(p) or p/3, where p is the number of predictors
                 max_depth=NULL, nodesize=5, max_nodes=NULL, # tree constraints
                 # maximum tree depth, minimum leaf size, and optional maximum number of leaf nodes
                 split_rule=c("SSE","Wilcoxon","pinball"), # SSE, rank-based Wilcoxon, or pinball loss splitting
                 replace = FALSE, # default is subsampling
                 subsample_prop= 1-exp(-1), # subsampling proportion; default is approximately 0.632
                 pinball_taus=c(0.1, 0.5, 0.9) # quantile levels used only for pinball splitting
) { 
  # Match and validate the requested splitting rule.
  split_rule <- match.arg(split_rule)

  # Basic argument checks. These do not change valid inputs, but give clearer
  # errors for values that would otherwise fail later inside sampling/tree growth.
  if (!is.numeric(ntree) || length(ntree) != 1L || is.na(ntree) ||
      !is.finite(ntree) || ntree < 1 || ntree != floor(ntree)) {
    stop("ntree must be a positive integer.")
  }
  
  if (!is.numeric(mtry) || length(mtry) != 1L || is.na(mtry) ||
      !is.finite(mtry) || mtry < 1 || mtry != floor(mtry)) {
    stop("mtry must be a positive integer.")
  }
  
  if (!is.numeric(nodesize) || length(nodesize) != 1L || is.na(nodesize) ||
      !is.finite(nodesize) || nodesize < 1 || nodesize != floor(nodesize)) {
    stop("nodesize must be a positive integer.")
  }
  
  if (!is.null(max_depth) &
      (!is.numeric(max_depth) || length(max_depth) != 1L || is.na(max_depth) ||
       !is.finite(max_depth) || max_depth < 1 || max_depth != floor(max_depth))) {
    stop("max_depth must be NULL or a positive integer.")
  }
  
  if (!is.null(max_nodes) &
      (!is.numeric(max_nodes) || length(max_nodes) != 1L || is.na(max_nodes) ||
       !is.finite(max_nodes) || max_nodes < 1 || max_nodes != floor(max_nodes))) {
    stop("max_nodes must be NULL or a positive integer.")
  }
  
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("replace must be TRUE or FALSE.")
  }
  
  if (!replace &
      (!is.numeric(subsample_prop) || length(subsample_prop) != 1L ||
       is.na(subsample_prop) || !is.finite(subsample_prop) ||
       subsample_prop <= 0 || subsample_prop > 1)) {
    stop("subsample_prop must be in (0, 1] when replace = FALSE.")
  }

  mf <- model.frame(formula, data, na.action = na.pass)
  yname <- names(mf)[1]
  xnames <- names(mf)[-1]
  if (!length(xnames)) stop("No predictors found.")
  
  if (anyNA(mf)) {
    stop("NAs detected. Please impute/drop before training.") # this implementation does not impute missing values
  }
  
  y <- as.numeric(mf[[1]])
  X <- mf[-1]
  n <- nrow(X)
  if (n < 1L) stop("Training data contain no observations.")
  
  if (split_rule == "SSE" || split_rule == "pinball") { # SSE and pinball splitting require finite numeric values.
    num_cols <- vapply(mf, is.numeric, logical(1))
    if (any(num_cols) & any(sapply(mf[num_cols], function(z) any(is.infinite(z))))) {
      stop("Inf/-Inf detected. Please impute/drop before training.")
    }
  }
  
  if (split_rule == "Wilcoxon") { # Rank-based splitting can rank Inf/-Inf values, but they often indicate data problems.
    num_cols <- vapply(mf, is.numeric, logical(1))
    if (any(num_cols) & any(sapply(mf[num_cols], function(z) any(is.infinite(z))))) {
      warning("Inf/-Inf detected. Although rank split can handle Infs and will proceed, please consider imputation/dropping before training.")
    }
  }
  
  y_split <- y
  
  if (split_rule == "pinball") {
    if (!is.numeric(pinball_taus) || length(pinball_taus) == 0L ||
        anyNA(pinball_taus) || any(!is.finite(pinball_taus))) {
      stop("pinball_taus must contain finite numeric values.")
    }
    pinball_taus <- sort(unique(pinball_taus))
    if (any(pinball_taus <= 0 | pinball_taus >= 1)) { # quantile levels must be valid probabilities
      stop("pinball_taus must be strictly between 0 and 1.")
    }
  }
  
  # Bootstrap if replace = TRUE; otherwise draw a subsample without replacement.
  nboot <- ifelse(replace, n, ceiling(subsample_prop*n))
  
  RF <- vector("list", ntree)
  
  # Store predictor types once so the recursive splitter can reuse them.
  x_is_num <- vapply(X, is.numeric, logical(1))
  
  for (t in seq_len(ntree)) {
    if (replace) { # bootstrap
      inbag_idx <- sample.int(n, n, replace=TRUE)
    } else { # subsampling
      inbag_idx <- sample.int(n, nboot, replace = FALSE)
    }
    
    # We pass 'y_split' for growing logic, but the tree indices link back to original 'y'
    RF[[t]] <- myRF_grow_tree(
      y = y_split, X = X, x_is_num = x_is_num,
      mtry = mtry, nodesize = nodesize, 
      split_rule = split_rule, pinball_taus = pinball_taus,
      max_depth = max_depth, max_nodes = max_nodes, 
      inbag_idx = inbag_idx)
  }
  
  # Return object stores original y for prediction
  structure(list(RF=RF, y=y, X=X, yname=names(mf)[1], xnames=xnames,
                 ntree=ntree, nodesize=nodesize, split_rule=split_rule,
                 pinball_taus=pinball_taus),
            class="myRF") # return as an "myRF" project
  
}


# Grow one recursive binary tree on the selected in-bag observations.
myRF_grow_tree <- function(y, X, x_is_num, mtry, nodesize, split_rule, pinball_taus,
                            max_depth, max_nodes, inbag_idx) {

  # Store nodes in a flat list. Internal nodes point to their left/right child IDs.
  # A small environment is used here so the nested helper functions can update the
  # growing node list
  tree_state <- new.env(parent = emptyenv())
  tree_state$nodes <- list()
  tree_state$splits_done <- 0L

  alloc <- function(node) {
    node_id <- length(tree_state$nodes) + 1L
    tree_state$nodes[[node_id]] <- node
    node_id
  }

  if (is.null(max_nodes)) {
    max_splits <- Inf
  } else {
    max_splits <- max(0L, as.integer(max_nodes) - 1L)
  }

  p <- ncol(X)

  build_node <- function(idx, depth) {
    # Stop if the node is too small, the tree is too deep, or the split budget is used.
    if (length(idx) < 2*nodesize || depth >= max_depth || tree_state$splits_done >= max_splits) {
      return(alloc(list(is_leaf = TRUE, train_idx = idx)))
    }

    # Randomly choose candidate predictors at this node.
    covars <- sample.int(p, size = min(mtry, p), replace = FALSE)
    # Track the best split across the candidate predictors.
    best <- list(score = -Inf)
    yj <- y[idx]

    for (j in covars) {
      xj <- X[[j]][idx]

      # Search for the best split for this predictor under the selected rule.
      if (split_rule == "pinball") {
        sp <- .find_best_pinball_split(yj, xj, nodesize, pinball_taus, x_is_num[j])
      } else if (split_rule == "Wilcoxon") {
        sp <- .find_best_Wilcoxon_split(yj, xj, nodesize, x_is_num[j])
      } else {
        sp <- .find_best_SSE_split(yj, xj, nodesize, x_is_num[j])
      }

      if (sp$score > best$score) { # update the optimal split score
        best <- sp
        best$covar <- j # covariate with the best split
      }
    }

    # If no valid split is found, this node becomes a leaf/terminal node.
    if (is.null(best$split_type)) {
      return(alloc(list(is_leaf = TRUE, train_idx = idx)))
    }

    # Send observations to the left or right child according to the selected split.
    xj_best <- X[[best$covar]][idx]
    if (best$split_type == "numeric") {
      left_mask <- xj_best <= best$split_point
    } else {
      left_mask <- as.character(xj_best) %in% best$split_point
    }

    # Count successful splits so max_nodes can be enforced.
    tree_state$splits_done <- tree_state$splits_done + 1L

    # Reserve a position for the current internal node. Children are built next,
    # and then the reserved position is filled with the split information.
    node_id <- alloc(list())
    left_id <- build_node(idx[left_mask], depth+1L)
    right_id <- build_node(idx[!left_mask], depth+1L)

    tree_state$nodes[[node_id]] <- list(
      is_leaf = FALSE, covar=best$covar, split_point=best$split_point,
      split_type=best$split_type, left = left_id, right = right_id)
    return(node_id)
  }

  build_node(inbag_idx, 0)
  return(list(nodes=tree_state$nodes, inbag=inbag_idx))
}


# Midpoint for a numeric split. For ordinary finite values this is identical
# to (a + b) / 2; the fallbacks only handle overflow or infinite endpoints.
.split_midpoint <- function(a, b) {
  mid <- (a + b) / 2
  if (is.finite(mid)) return(mid)

  if (is.finite(a) & is.finite(b)) {
    return(a / 2 + b / 2)
  }
  if (is.finite(a) & is.infinite(b)) {
    return(a)
  }
  if (is.infinite(a) & is.finite(b)) {
    return(a)
  }
  if (is.infinite(a) & is.infinite(b) & a < b) {
    return(a)
  }
  mid
}


# Find the best SSE split for one predictor within a node.
.find_best_SSE_split <- function(y, x, nodesize, is_num) {
  n <- length(y)
  if (n < 2*nodesize) return(list(score = -Inf))
  
  if (is_num) {
    ord <- order(x)
    xs <- x[ord]
    ys <- y[ord]
    
    # Identify valid split indices where x changes
    # Use 'diff' to find boundaries
    valid_cuts <- which(diff(xs) != 0)
    # Filter for required minimum node size
    valid_cuts <- valid_cuts[valid_cuts >= nodesize & valid_cuts <= (n - nodesize)]
    
    if (length(valid_cuts) == 0) return(list(score = -Inf))
    
    # Parent SSE = sum(y^2) - sum(y)^2/n (we only care about maximizing reduction)
    # Reduction = (sumL^2/nL + sumR^2/nR) - (sumTotal^2/n)
    # So essentially, we just maximize (sumL^2/nL + sumR^2/nR)
    cy <- cumsum(ys)
    total_y <- cy[n]
    sumL <- cy[valid_cuts]
    nL <- valid_cuts
    sumR <- total_y - sumL
    nR <- n - nL
    
    # Maximizing this quantity is equivalent to maximizing SSE reduction.
    crit <- (sumL^2 / nL) + (sumR^2 / nR)
    
    best_idx <- which.max(crit)
    k <- valid_cuts[best_idx]
    
    return(list(
      score = crit[best_idx],
      split_point = .split_midpoint(xs[k], xs[k+1]), # split at midpoint
      split_type = "numeric"
    ))
    
  } else { # for categorical covariate X
    # sort categories by mean(y)
    ux <- unique(x)
    if (length(ux) < 2) return(list(score = -Inf))
    
    # Aggregate outcome sums and counts by category.
    x_fac <- factor(x)
    levs <- levels(x_fac)
    sums <- as.vector(tapply(y, x_fac, sum))
    cnts <- as.vector(tapply(y, x_fac, length))
    
    # Filter empty levels
    keep <- cnts > 0
    sums <- sums[keep]; cnts <- cnts[keep]; levs <- levs[keep]
    
    if (length(levs) < 2) return(list(score = -Inf))
    
    means <- sums / cnts
    ord <- order(means)
    sums <- sums[ord]; cnts <- cnts[ord]; levs <- levs[ord]
    
    # Search binary splits after ordering categories by their node-specific mean outcome.
    cs_sum <- cumsum(sums)
    cs_cnt <- cumsum(cnts)
    
    total_y <- cs_sum[length(cs_sum)]
    n <- cs_cnt[length(cs_cnt)]
    
    # Candidates
    valid_idx <- which(cs_cnt >= nodesize & (n - cs_cnt) >= nodesize)
    if (length(valid_idx) == 0) return(list(score = -Inf))
    
    # Do not allow the split that sends all categories to the left child.
    valid_idx <- valid_idx[valid_idx < length(levs)]
    if (length(valid_idx) == 0) return(list(score = -Inf))
    
    sumL <- cs_sum[valid_idx]
    nL <- cs_cnt[valid_idx]
    sumR <- total_y - sumL
    nR <- n - nL
    
    crit <- (sumL^2 / nL) + (sumR^2 / nR)
    best_i <- which.max(crit)
    k <- valid_idx[best_i]
    
    return(list(score = crit[best_i],
                split_point = levs[1:k],
                split_type = "categorical"))
  }
}


# Find the best Wilcoxon split
.find_best_Wilcoxon_split <- function(y, x, nodesize, is_num) {
  n <- length(y)
  if (n < 2*nodesize) return(list(score = -Inf))
  
  # Calculate LOCAL midranks (for handling ties) for the current node
  r <- rank(y, ties.method = "average") # midrank
  
  # Pre-calculate constants for the statistic
  # Expected sum of ranks for a group of size k is k * (n + 1) / 2
  grand_mean_rank <- (n + 1) / 2
  
  if (is_num) {
    # Sort data by predictor X
    ord <- order(x)
    xs <- x[ord]
    rs <- r[ord] # Sorted local ranks
    
    valid_cuts <- which(diff(xs) != 0)
    valid_cuts <- valid_cuts[valid_cuts >= nodesize & valid_cuts <= (n - nodesize)]
    
    if (length(valid_cuts) == 0) return(list(score = -Inf))
    
    # cum sum of ranks
    cum_rank <- cumsum(rs)
    sumL <- cum_rank[valid_cuts]
    nL <- valid_cuts # Size of left node
    
    # Wilcoxon split score.  For a left child of size nL, the expected
    # rank sum under no separation is nL * (n + 1) / 2.  We maximize the
    # squared standardized deviation, similar to SSE split
    numerator <- (sumL - nL * grand_mean_rank)^2
    nR <- n - nL
    crit <- numerator / (nL * nR)
    
    best_idx <- which.max(crit)
    k <- valid_cuts[best_idx]
    
    return(list(score = crit[best_idx],
                split_point = .split_midpoint(xs[k], xs[k+1]),
                split_type = "numeric"))
    
  } else { # for categorical covariate X
    ux <- unique(x)
    if (length(ux) < 2) return(list(score = -Inf))
    
    x_fac <- factor(x)
    levs <- levels(x_fac)
    
    # For categorical, we need to order the categories to allow for finding the split.
    # order by mean(midrank(Y)).
    sum_ranks <- as.vector(tapply(r, x_fac, sum))
    cnts <- as.vector(tapply(r, x_fac, length))
    
    keep <- cnts > 0
    sum_ranks <- sum_ranks[keep]
    cnts <- cnts[keep]
    levs <- levs[keep]
    
    if (length(levs) < 2) return(list(score = -Inf))
    
    # Sort categories by level-wise average rank
    avg_ranks <- sum_ranks / cnts
    ord <- order(avg_ranks)
    
    sum_ranks <- sum_ranks[ord]
    cnts <- cnts[ord]
    levs <- levs[ord]
    
    # Search splits along sorted category levels
    cs_sum_rank <- cumsum(sum_ranks)
    cs_cnt <- cumsum(cnts)
    
    valid_idx <- which(cs_cnt >= nodesize & (n - cs_cnt) >= nodesize)
    valid_idx <- valid_idx[valid_idx < length(levs)]
    
    if (length(valid_idx) == 0) return(list(score = -Inf))
    
    sumL <- cs_sum_rank[valid_idx]
    nL <- cs_cnt[valid_idx]
    nR <- n - nL
    
    numerator <- (sumL - nL * grand_mean_rank)^2
    crit <- numerator / (nL * nR)
    
    best_i <- which.max(crit)
    k <- valid_idx[best_i]
    
    return(list(
      score = crit[best_i],
      split_point = levs[1:k],
      split_type = "categorical"))
  }
}


# Helpers for pinball-loss-based splitting.
.pinball_loss_vec <- function(e, tau) (tau - (e < 0)) * e

.pinball_node_loss <- function(y, taus) {
  # The empirical quantile calculation below assumes sorted outcomes.
  # Sorting here is safer because callers may pass outcomes in predictor order.
  y <- sort(y)
  s <- 0
  n <- length(y)
  
  for(tau in taus) {
    # Type-7 empirical quantile, same as R's default quantile() convention.
    h <- (n - 1) * tau + 1
    h_floor <- floor(h)
    
    if (h_floor == h) {
      q <- y[h_floor]
    } else {
      # Linear interpolation between the two nearest indices
      q <- y[h_floor] + (h - h_floor) * (y[h_floor + 1] - y[h_floor])
    }
    
    s <- s + sum(.pinball_loss_vec(y - q, tau))
  }
  return(s)
}


# Find the best approximate pinball-loss split for one predictor within a node.
# Numeric predictors use a grid of candidate split points for speed; categorical
# predictors are ordered by median outcome and then scanned.
.find_best_pinball_split <- function(y, x, nodesize, taus, is_num, max_cuts = 20) {
  n <- length(y)
  if (n < 2 * nodesize) return(list(score = -Inf))
  
  # Use a coarser candidate grid in large nodes to reduce computation time.
  if (n > 200) max_cuts <- 40
  
  parent_pbloss <- .pinball_node_loss(y, taus)
  best <- list(score = -Inf)
  
  if (is_num) {
    ord <- order(x)
    xs <- x[ord]; ys <- y[ord]
    
    # Candidate split positions must leave at least nodesize observations on each side.
    valid_range <- nodesize:(n - nodesize)
    if (length(valid_range) == 0) return(best)
    
    # Evaluate either all valid positions or an evenly spaced subset.
    if (length(valid_range) > max_cuts) {
      cand_idx <- unique(round(seq(min(valid_range), max(valid_range), length.out = max_cuts)))
    } else {
      cand_idx <- valid_range
    }
    # A numeric split is valid only where adjacent ordered predictor values differ.
    cand_idx <- cand_idx[xs[cand_idx] != xs[cand_idx+1]]
    
    if (length(cand_idx) == 0) return(best)
    
    for (k in cand_idx) {
      pbloss <- .pinball_node_loss(ys[1:k], taus) + .pinball_node_loss(ys[(k+1):n], taus)
      reduction <- parent_pbloss - pbloss
      if (reduction > best$score) {
        best <- list(score=reduction, split_point=.split_midpoint(xs[k], xs[k+1]), split_type="numeric")
      }
    }
  } else {
    # For categorical predictors, order categories by their median outcome and scan splits.
    x_fac <- factor(x)
    if (length(levels(x_fac)) < 2) return(best)
    
    # Order categories by the median outcome as a practical heuristic.
    lev_meds <- tapply(y, x_fac, median)
    lev_cnts <- tapply(y, x_fac, length)
    valid_levs <- names(lev_cnts)[lev_cnts > 0]
    if (length(valid_levs) < 2) return(best)
    
    ord <- order(lev_meds[valid_levs])
    sorted_levs <- valid_levs[ord]
    
    K <- length(sorted_levs)
    # Scan binary splits along the ordered categories.
    for (k in 1:(K-1)) {
      left_set <- sorted_levs[1:k]
      mask <- x_fac %in% left_set
      nL <- sum(mask); nR <- n - nL
      if (nL < nodesize || nR < nodesize) next
      
      pbloss <- .pinball_node_loss(y[mask], taus) + .pinball_node_loss(y[!mask], taus)
      reduction <- parent_pbloss - pbloss
      if (reduction > best$score) {
        best <- list(score=reduction, split_point=left_set, split_type="categorical")
      }
    }
  }
  return(best)
}


# Compute the forest weight matrix for new observations.
# Row i contains the weights assigned to training outcomes for prediction at newdata[i, ].
# The subsample/bootstrap sample is used to construct each tree. Once the tree is fixed,
# all original training observations are passed through the tree and observations in the
# same terminal leaf contribute once to that tree's empirical distribution.
.qrf_get_weights <- function(model, newdata) {
  ntree <- model$ntree
  ntrain <- length(model$y)
  n_new <- nrow(newdata)

  # Store sparse-matrix coordinates for each tree, then combine them at the end.
  i_list <- vector("list", ntree)
  j_list <- vector("list", ntree)
  x_list <- vector("list", ntree)

  for (t in seq_len(ntree)) {
    tree <- model$RF[[t]]

    # Route all new observations through the fitted tree.
    leaf_ids_new <- .get_leaf_ids_vectorized(tree, newdata)

    # After the partition is fixed, route all original training observations
    # through the tree. Each training observation contributes once to the
    # terminal-node empirical distribution, regardless of whether it was in-bag.
    train_ids <- if (!is.null(tree$leaf_id_by_obs)) {
      tree$leaf_id_by_obs
    } else {
      .get_leaf_ids_vectorized(tree, model$X)
    }

    new_by_leaf <- split(seq_len(n_new), leaf_ids_new)
    train_by_leaf <- split(seq_len(ntrain), train_ids)

    common_leaves <- intersect(names(new_by_leaf), names(train_by_leaf))

    # Temporary lists for this specific tree
    i_leaf <- vector("list", length(common_leaves))
    j_leaf <- vector("list", length(common_leaves))
    x_leaf <- vector("list", length(common_leaves))

    for (k in seq_along(common_leaves)) {
      lid <- common_leaves[k]
      ii_new <- new_by_leaf[[lid]]
      ii_train <- train_by_leaf[[lid]]

      nL <- length(ii_train)
      if (nL > 0) {
        # Each new observation receives equal weight on all original training
        # observations in the same fitted leaf.
        n_combos <- length(ii_new) * nL
        i_leaf[[k]] <- rep(ii_new, each = nL)
        j_leaf[[k]] <- rep(ii_train, times = length(ii_new))
        x_leaf[[k]] <- rep(1 / nL, n_combos)
      }
    }

    # Store this tree's sparse-matrix entries.
    i_list[[t]] <- unlist(i_leaf, use.names = FALSE)
    j_list[[t]] <- unlist(j_leaf, use.names = FALSE)
    x_list[[t]] <- unlist(x_leaf, use.names = FALSE)
  }

  # Combine entries across trees.
  i_all <- unlist(i_list, use.names = FALSE)
  j_all <- unlist(j_list, use.names = FALSE)
  x_all <- unlist(x_list, use.names = FALSE)

  # This is where Matrix package comes in!!
  # Build a sparse matrix; duplicate entries are summed automatically.
  final_weight <- Matrix::sparseMatrix(
    i = i_all, # rows
    j = j_all, # columns
    x = x_all,
    dims = c(n_new, ntrain))

  # Average the tree-level weights over the forest.
  final_weight <- final_weight / ntree
  return(final_weight)
}

.qrf_quantile_helper <- function(wi, ysorted, probs) {
  weight_sum <- sum(wi)
  if (weight_sum <= 0) return(rep(NA_real_, length(probs)))
  
  w_norm <- wi / weight_sum
  ecdff <- cumsum(w_norm)
  
  sapply(probs, function(p) {
    if (p <= 0) return(ysorted[1L])
    if (p >= 1) return(ysorted[length(ysorted)])
    
    # Binary search using findInterval (left.open=TRUE mimics >= behavior)
    idx <- findInterval(p, ecdff, left.open = TRUE) + 1L
    
    # Fallback to prevent indexing out of bounds on floating point errors
    if (idx > length(ysorted)) idx <- length(ysorted)
    
    return(ysorted[idx])
  })
}


# Step-CDF probabilities from WEIGHTS (classic QRF): vectorized over outcomes
.qrf_prob_from_weights <- function(wi, y_sorted, outcomes, lower.tail = TRUE, include = TRUE) {
  s <- sum(wi)
  if (s <= 0) return(rep(NA_real_, length(outcomes)))
  w <- wi / s
  Fw <- cumsum(w)  # F(y_k) = P(Y <= y_k)
  
  # counts of grid points strictly below (<) and at-or-below (<=) each outcome
  idx_le <- findInterval(outcomes, y_sorted) # {y_i <= y}
  idx_lt <- findInterval(outcomes, y_sorted, left.open = TRUE) # {y_i <  y}
  
  F_le <- ifelse(idx_le > 0L, Fw[idx_le], 0)
  F_lt <- ifelse(idx_lt > 0L, Fw[idx_lt], 0)
  
  if (lower.tail) {
    if (include) F_le else F_lt
  } else {
    if (include) 1 - F_lt else 1 - F_le
  }
}


# Weighted quantiles with optional linear interpolation between support points.
.qrf_quantile_linear_helper <- function(wi, ysorted, probs) {
  weight_sum <- sum(wi)
  if (weight_sum <= 0) return(rep(NA_real_, length(probs)))

  w_norm <- wi / weight_sum
  ecdff <- cumsum(w_norm)

  # A degenerate weighted distribution has only one support point with positive mass.
  # approx() needs at least two distinct x-values, so handle that case directly.
  positive_idx <- which(w_norm > 0)
  if (length(positive_idx) == 1L) {
    return(rep(ysorted[positive_idx], length(probs)))
  }

  # Use approx() to find the y value for a given probability (x=CDF, y=Value), kinda like linear interpolation
  # rule=2 ensures we return min/max y if prob is outside observed CDF range (0 or 1)
  # ties="ordered" assumes inputs are already sorted for efficiency
  approx(x = ecdff, y = ysorted, xout = probs, method = "linear", rule = 2, ties = "ordered")$y
}

# Interpolated probabilities from a weighted empirical CDF.
.qrf_prob_linear_from_weights <- function(wi, ysorted, outcomes, lower.tail = TRUE) {
  weight_sum <- sum(wi)
  if (weight_sum <= 0) return(rep(NA_real_, length(outcomes)))

  w_norm <- wi / weight_sum
  ecdff <- cumsum(w_norm)

  # For interpolation, we need unique x-coordinates (unique y values).
  # We take the cumulative probability at the *last* occurrence of each unique y
  # (i.e., the full CDF value at that point).
  keep_idx <- !duplicated(ysorted, fromLast = TRUE)
  y_unique <- ysorted[keep_idx]
  f_unique <- ecdff[keep_idx]

  if (length(y_unique) == 1L) {
    probs <- as.numeric(outcomes >= y_unique[1L])
    if (!lower.tail) return(1 - probs)
    return(probs)
  }

  # Interpolate F(y) (x=Value, y=CDF)
  probs <- approx(x = y_unique, y = f_unique, xout = outcomes, 
                  method = "linear", rule = 2, ties = "ordered")$y

  # Correction: approx with rule=2 extends the min F value (e.g. F(y_1)) to inputs < y_1
  # For probabilities, inputs < min(y) should be 0.
  probs[outcomes < y_unique[1L]] <- 0

  if (!lower.tail) {
    return(1 - probs)
  }
  return(probs)
}

# Helper function to calculate continuous ranked probability score (CRPS)
.calc_crps <- function(W, y_train, y_obsvd) {
  # W: Sparse Weight matrix (n_test x n_train)
  # y_train: Training outcomes (length n_train)
  # y_obsvd: Observed outcomes for the test set (length n_test)
  
  if (length(y_obsvd) != nrow(W)) {
    stop("Length of 'y_obsvd' must match the number of rows in 'W'.")
  }
  n_test <- nrow(W)
  n_train <- length(y_train)
  
  # Sort weights by Y
  ord <- order(y_train)
  y_sorted <- y_train[ord]
  W_sorted <- W[, ord, drop = FALSE]
  # collpase weights for unique Y values
  is_unique <- !duplicated(y_sorted)
  y_unique <- y_sorted[is_unique]
  K <- length(y_unique)
  
  # Use the sparse indicator matrix trick for aggregation
  if (K < n_train) {
    # Create group mapping (1 to K)
    group_ids <- cumsum(is_unique)
    
    # Build a sparse mapping matrix (n_train x K)
    # Again, use Matrix::sparseMatrix()
    indicator_matrix <- Matrix::sparseMatrix(
      i = seq_len(n_train),
      j = group_ids,
      x = 1,
      dims = c(n_train, K))
    
    # Fast sparse-sparse matrix multiplication to aggregate columns!
    W_agg <- W_sorted %*% indicator_matrix
  } else {
    W_agg <- W_sorted
  }
  
  # Pre-calculate differences between unique support points for Term 2
  dy <- diff(y_unique)
  
  # Output vector
  crps_out <- numeric(n_test)
  
  # 3. Row-by-Row CRPS calculation to avoid dense matrices
  for (i in seq_len(n_test)) {
    
    # Extract row i as a standard numeric vector
    wi <- as.numeric(W_agg[i, ])
    
    # Other than the integrated Brier scores indicated in the main text,
    # CRPS can also be expressed as E|Z - y_obs| - 0.5 * E|Z - Z'|,
    # where Z, Z' are iid RVs following F. For convenience, we use this expression instead.
    # Term 1: E|Z - y_obs|
    # Vectorized distance of support points from the single observation
    dist_i <- abs(y_unique - y_obsvd[i])
    term1 <- sum(wi * dist_i)
    
    # Term 2: 0.5 * E|Z - Z'| = integral F(z)(1-F(z)) dz
    if (K > 1) {
      Fi <- cumsum(wi)
      Fi[K] <- 1.0 # in case of floating point precision errors
      
      term2 <- sum(Fi[-K] * (1 - Fi[-K]) * dy)
    } else {
      term2 <- 0
    }
    
    # CRPS = Term 1 - Term 2
    crps_out[i] <- term1 - term2
  }
  
  return(crps_out)
}


# Calculate the Cramer-von Mises statistic from PIT values, or return probability-scale residual (PSR) values.
.calc_pit_cvm <- function(W, y_train, y_obsvd, PSR) {
  if (length(y_obsvd) != nrow(W)) {
    stop("Length of 'y_obsvd' must match the number of rows in 'W'.")
  }
  
  n_test <- nrow(W)
  n_train <- length(y_train)
  
  # Sort weights by unique Y
  ord <- order(y_train)
  y_sorted <- y_train[ord]
  W_sorted <- W[, ord, drop = FALSE]
  
  is_unique <- !duplicated(y_sorted)
  y_unique <- y_sorted[is_unique]
  K <- length(y_unique)
  
  # That sparse indicator matrix trick
  if (K < n_train) {
    group_ids <- cumsum(is_unique)
    
    # Build a sparse mapping matrix to cleanly aggregate duplicate Y values
    indicator_matrix <- Matrix::sparseMatrix(
      i = seq_len(n_train),
      j = group_ids,
      x = 1,
      dims = c(n_train, K))
    
    # Fast sparse matrix multiplication to sum groups natively
    W_agg <- W_sorted %*% indicator_matrix
  } else {
    W_agg <- W_sorted
  }
  
  # Pre-calculate the outcome intervals for the whole test set
  idx_le <- findInterval(y_obsvd, y_unique)  # P(Y <= y_obs)
  idx_lt <- findInterval(y_obsvd, y_unique, left.open = TRUE) # P(Y < y_obs)
  
  U_vec <- numeric(n_test)
  
  # Row-by-row PIT calculation (avoid building a massive F_agg matrix)
  for (i in seq_len(n_test)) {
    
    # Extract the sparse row as a standard numeric vector
    wi <- as.numeric(W_agg[i, ])
    
    # Fast cumulative sum for this single observation
    Fi <- cumsum(wi)
    
    # Extract P(Y <= y_obs) and P(Y < y_obs)
    p_le_i <- if (idx_le[i] > 0) Fi[idx_le[i]] else 0
    p_lt_i <- if (idx_lt[i] > 0) Fi[idx_lt[i]] else 0
    
    # P(Y = y_obs) is the difference
    p_eq_i <- p_le_i - p_lt_i
    
    # PIT value for observation i
    U_vec[i] <- p_lt_i + 0.5 * p_eq_i
  }
  
  # Return PSR or CvM statistic
  if (PSR) {
    return(2 * U_vec - 1) # rescale to PSR
  } else {
    # Calculate CvM Statistic
    U_sorted <- sort(U_vec)
    seq_i <- seq_len(n_test)
    term_sum <- sum((U_sorted - (2 * seq_i - 1) / (2 * n_test))^2)
    W2 <- term_sum + 1 / (12 * n_test)
    
    return(W2)
  }
}


# Probability-scale residuals (PSRs) for ordinal discrete (or mixed discrete-continuous) outcomes
# Based on Shepherd, Li & Liu (2016), Section 4:
#   PSR = 2*F*(y) - f*(y) - 1 = F*(y-) + F*(y) - 1 = P(Y* < y) - P(Y* > y)
# Properties under correct specification:
#   E(PSR) = 0;  Var(PSR) = {1 - sum f*(y_k)^3} / 3
.calc_psr_discrete <- function(W, y_train, y_obsvd) {
  if (length(y_obsvd) != nrow(W)) {
    stop("Length of 'y_obsvd' must match the number of rows in 'W'.")
  }
  
  n_obs <- nrow(W)
  n_train <- length(y_train)
  
  # Sort and aggregate weights to unique Y values
  ord <- order(y_train)
  y_sorted <- y_train[ord]
  W_sorted <- W[, ord, drop = FALSE]
  
  is_unique <- !duplicated(y_sorted)
  y_unique <- y_sorted[is_unique]
  K <- length(y_unique)
  
  if (K < n_train) {
    group_ids <- cumsum(is_unique)
    ind_mat <- Matrix::sparseMatrix(
      i = seq_len(n_train), j = group_ids, x = 1,
      dims = c(n_train, K))
    W_agg <- W_sorted %*% ind_mat
  } else {
    W_agg <- W_sorted
  }
  
  # Precompute interval indices
  idx_le <- findInterval(y_obsvd, y_unique) # F*(y)
  idx_lt <- findInterval(y_obsvd, y_unique, left.open = TRUE) # F*(y-)
  
  F_le <- numeric(n_obs)
  F_lt <- numeric(n_obs)
  
  for (i in seq_len(n_obs)) {
    wi <- as.numeric(W_agg[i, ])
    Fi <- cumsum(wi)
    F_le[i] <- if (idx_le[i] > 0L) Fi[idx_le[i]] else 0
    F_lt[i] <- if (idx_lt[i] > 0L) Fi[idx_lt[i]] else 0
  }
  
  f_eq <- F_le - F_lt # f*(y) = P(Y* = y)
  psr <- F_lt + F_le - 1 # PSR = F*(y-) + F*(y) - 1
  
  return(list(psr = psr, F_le = F_le, F_lt = F_lt, f_eq = f_eq))
}


# Predict from a fitted myRF object.
#
# what = "quantiles" returns weighted empirical quantiles.
# what = "probability" returns threshold probabilities at the values in outcomes.
# what = "CRPS", "CvM", "PSR", and "PSR discrete" require observed outcomes
# through y_obsvd and are mainly used for evaluation.
predict.myRF <- function(object, newdata,
                         probs=c(0.1, 0.5, 0.9),
                         what=c("quantiles", "mean", "median", "ECDF", "weights", "probability", "CRPS", "CvM", "PSR", "PSR discrete"),
                         outcomes=NULL,
                         y_obsvd=NULL, 
                         lower.tail = TRUE, include = TRUE, 
                         linear_intpl = FALSE, ...) {
  
  what <- match.arg(what)
  
  newdata <- as.data.frame(newdata)[, object$xnames, drop=FALSE]
  if (anyNA(newdata)) stop("newdata contains NA(s) at prediction time.")
  
  # Calculate forest weights
  W <- .qrf_get_weights(object, newdata)
  
  y <- object$y
  ord <- order(y)
  ysorted <- y[ord]
  Wreord <- W[, ord, drop=FALSE]
  
  if (what == "CRPS") {
    if (is.null(y_obsvd)) stop("You must provide 'y_obsvd' (observed y) for CRPS.")
    return(.calc_crps(W, y, y_obsvd))
  }
  
  if (what == "CvM" || what == "PSR") {
    if (is.null(y_obsvd)) stop("You must provide 'y_obsvd' (observed y) for CvM.")
    if (what == "CvM") {
      return(.calc_pit_cvm(W, y, y_obsvd, PSR = FALSE))
    } else {
      return(.calc_pit_cvm(W, y, y_obsvd, PSR = TRUE))
    }
  }
  
  if (what == "PSR discrete") {
    if (is.null(y_obsvd)) stop("You must provide 'y_obsvd' (observed y) for PSR discrete.")
    return(.calc_psr_discrete(W, y, y_obsvd))
  }
  
  if (what == "mean") {
    return(as.numeric(Wreord %*% ysorted))
  }
  
  if (what == "median") {
    if (linear_intpl) {
      qs <- apply(Wreord, 1, .qrf_quantile_linear_helper, ysorted=ysorted, probs=0.5)
    } else {
      qs <- apply(Wreord, 1, .qrf_quantile_helper, ysorted=ysorted, probs=0.5)
    }
    return(as.numeric(qs))
  }
  
  if (what == "quantiles") {
    if (linear_intpl) {
      Q <- t(apply(Wreord, 1, .qrf_quantile_linear_helper, ysorted=ysorted, probs=probs))
    } else {
      Q <- t(apply(Wreord, 1, .qrf_quantile_helper, ysorted=ysorted, probs=probs))
    }
    
    if (length(probs)==1) {
      out = data.frame(matrix(data=Q[1,],ncol=1))
    } else {
      out <- as.data.frame(Q)
    }
    colnames(out) <- paste0("q_", probs)
    return(out)
  }
  
  if (what == "ECDF" || what == "weights") {
    if (nrow(newdata) != 1) {
      stop("what='ECDF' and what='weights' only support nrow(newdata) == 1")
    }
    wi <- Wreord[1, ]
    weight_sum <- sum(wi)
    
    if (weight_sum > 0) {
      w_norm <- wi / weight_sum
    } else {
      w_norm <- wi 
    }
    
    if (what == "weights") {
      return(data.frame(y.sorted = ysorted, weight = w_norm))
    }
    if (what == "ECDF") {
      ecdff <- cumsum(w_norm)
      out <- data.frame(y.sorted = ysorted, ECDF = ecdff)
      if (!linear_intpl) {
        # Keep the step-function representation compact. For tied outcome values,
        # retain the last row so the full probability mass at that value is included.
        out <- out[!duplicated(out$y.sorted, fromLast = TRUE), , drop = FALSE]
        out <- out[!duplicated(out$ECDF), , drop = FALSE]
      }
      return(out)
    }
  }
  
  if (what == "probability") {
    if (is.null(outcomes)) stop("Provide 'outcomes' when what='probability'.")
    if (linear_intpl) {
      P0 <- apply(Wreord, 1, function(wi)
        .qrf_prob_linear_from_weights(wi, ysorted, outcomes, lower.tail = lower.tail))
    } else {
      P0 <- apply(Wreord, 1, function(wi)
        .qrf_prob_from_weights(wi, ysorted, outcomes,
                               lower.tail = lower.tail, include = include))
    }
    if (length(outcomes) == 1L) {
      P <- matrix(P0, ncol = 1L)
    } else {
      P <- t(P0)
    }
    P <- as.data.frame(P)
    colnames(P) <- as.character(outcomes)
    return(P)
  }
}


# Out-of-bag prediction from a fitted myRF object.
# For each subject, only trees in which that subject was out-of-bag contribute to the prediction.
# Under the current estimator, when response_members is omitted the terminal-node
# distribution uses all training observations other than the subject being predicted.
# response_members = "inbag" is retained as an explicit strict in-bag alternative.
predict_oob.myRF <- function(object,
                             probs = c(.1, .5, .9),
                             what = c("quantiles", "mean", "median", "probability", "CRPS", "CvM", "PSR", "PSR discrete"),
                             outcomes = NULL,
                             lower.tail = TRUE,
                             include = TRUE,
                             linear_intpl = FALSE,
                             response_members = c("inbag", "all"), ...) {
  
  what <- match.arg(what)
  # Keep the public function signature unchanged. For the current estimator, an
  # omitted response_members argument uses all training observations except the
  # OOB subject itself. Users may still request the historical strict in-bag option.
  if (missing(response_members)) {
    response_members <- "all"
  } else {
    response_members <- match.arg(response_members)
  }
  y <- object$y
  
  # OOB weights are built from the fitted forest. Under the current estimator,
  # each OOB subject is predicted from all other training observations in the
  # terminal node. This avoids using the subject's own outcome in its prediction.
  ow <- .oob_weights_fast(object, response_members = response_members)
  W <- ow$W
  denom <- ow$denom
  
  ok <- denom > 0L
  inv_denom <- rep(0, length(denom))
  inv_denom[ok] <- 1 / denom[ok]
  W <- Matrix::Diagonal(x = inv_denom) %*% W
  
  # CRPS or CvM or PSR
  if (what == "CRPS") {
    crps_out <- rep(NA_real_, length(y))
    if (any(ok)) {
      crps_out[ok] <- .calc_crps(as(W[ok, , drop=FALSE], "RsparseMatrix"), y, y[ok])
    }
    return(crps_out)
  }
  
  if (what == "CvM") {
    if (any(ok)) {
      return(.calc_pit_cvm(as(W[ok, , drop=FALSE], "RsparseMatrix"),
                           y, y[ok], PSR = FALSE))
    }
    return(NA_real_)
  }

  if (what == "PSR") {
    psr_out <- rep(NA_real_, length(y))
    if (any(ok)) {
      psr_out[ok] <- .calc_pit_cvm(as(W[ok, , drop=FALSE], "RsparseMatrix"),
                                   y, y[ok], PSR = TRUE)
    }
    return(psr_out)
  }
  
  if (what == "PSR discrete") {
    if (any(ok)) {
      result <- .calc_psr_discrete(as(W[ok, , drop=FALSE], "RsparseMatrix"), y, y[ok])
      # Expand back to full length with NAs for in-bag-only observations
      psr_full <- rep(NA_real_, length(y))
      F_le_full <- rep(NA_real_, length(y))
      F_lt_full <- rep(NA_real_, length(y))
      f_eq_full <- rep(NA_real_, length(y))
      psr_full[ok] <- result$psr
      F_le_full[ok] <- result$F_le
      F_lt_full[ok] <- result$F_lt
      f_eq_full[ok] <- result$f_eq
      return(list(psr = psr_full, F_le = F_le_full, F_lt = F_lt_full, f_eq = f_eq_full))
    }
    return(list(psr = rep(NA_real_, length(y)), F_le = rep(NA_real_, length(y)),
                F_lt = rep(NA_real_, length(y)), f_eq = rep(NA_real_, length(y))))
  }
  
  ord <- order(y)
  ys <- y[ord]
  Wre <- W[, ord, drop = FALSE]
  
  if (what == "mean") {
    mu <- as.numeric(Wre %*% ys)
    mu[!ok] <- NA_real_
    return(mu)
  }
  
  if (what %in% c("median", "quantiles")) {
    p_targets <- if (what == "median") 0.5 else probs
    out <- matrix(NA_real_, nrow = length(y), ncol = length(p_targets))
    
    if (any(ok)) {
      valid_idx <- which(ok)
      Wre_row <- as(Wre, "RsparseMatrix") # convert to RsparseMatrix object
      for (i in valid_idx) {
        wi <- as.numeric(Wre_row[i, ])
        if (linear_intpl) {
          out[i, ] <- .qrf_quantile_linear_helper(wi, ys, p_targets)
        } else {
          out[i, ] <- .qrf_quantile_helper(wi, ys, p_targets)
        }
      }
    }
    
    if (what == "median") return(as.vector(out))
    colnames(out) <- as.character(probs)
    return(as.data.frame(out))
  }
  
  if (what == "probability") {
    if (is.null(outcomes)) stop("Provide 'outcomes' when what='probability'.")
    P <- matrix(NA_real_, nrow = length(y), ncol = length(outcomes))
    
    if (any(ok)) {
      if (linear_intpl) {
        W_ok_row <- as(Wre[ok, , drop = FALSE], "RsparseMatrix")
        P[ok, ] <- t(apply(W_ok_row, 1, function(wi)
          .qrf_prob_linear_from_weights(wi, ys, outcomes, lower.tail = lower.tail)))
      } else {
        # Fully vectorized step-CDF evaluation
        idx_le <- findInterval(outcomes, ys)
        idx_lt <- findInterval(outcomes, ys, left.open = TRUE)
        
        W_ok <- Wre[ok, , drop = FALSE]
        F_ok <- W_ok
        for (j in 2:ncol(F_ok)) {
          F_ok[, j] <- F_ok[, j-1] + F_ok[, j]
        }
        
        for (k in seq_along(outcomes)) {
          if (lower.tail) {
            idx <- if (include) idx_le[k] else idx_lt[k]
            P[ok, k] <- if (idx > 0) F_ok[, idx] else 0
          } else {
            idx <- if (include) idx_lt[k] else idx_le[k]
            P[ok, k] <- if (idx > 0) 1 - F_ok[, idx] else 1
          }
        }
      }
    }
    colnames(P) <- as.character(outcomes)
    return(as.data.frame(P))
  }
}

# Helper function to obtain CRPS using OOB weights
.oob_crps_from_weights <- function(object, ow) {
  y <- object$y
  W <- ow$W
  denom <- ow$denom
  ok <- denom > 0L
  
  crps_out <- rep(NA_real_, length(y))
  if (!any(ok)) return(crps_out)
  
  inv_denom <- rep(0, length(denom))
  inv_denom[ok] <- 1 / denom[ok]
  W <- Matrix::Diagonal(x = inv_denom) %*% W
  crps_out[ok] <- .calc_crps(as(W[ok, , drop = FALSE], "RsparseMatrix"), y, y[ok])
  crps_out
}


# Build the OOB weight matrix for one fixed tree. This helper is used only for
# tree-level permutation importance, where CRPS is calculated separately within
# each tree before the losses are averaged across trees.
.tree_oob_weights_for_vimp <- function(object, tree_index, response_members,
                                       lids_base, used_vars,
                                       permute_var = NULL) {
  n <- length(object$y)
  X <- object$X
  tree <- object$RF[[tree_index]]
  inbag <- tree$inbag

  is_oob <- rep(TRUE, n)
  is_oob[unique(inbag)] <- FALSE
  oob_idx <- which(is_oob)
  if (length(oob_idx) == 0L) {
    return(list(W = NULL, oob_idx = oob_idx))
  }

  if (is.null(permute_var) || !(permute_var %in% used_vars)) {
    oob_lids <- lids_base[oob_idx]
  } else {
    x_perm <- X[[permute_var]]
    x_perm[oob_idx] <- sample(x_perm[oob_idx], length(oob_idx), replace = FALSE)
    oob_lids <- .get_leaf_ids_vectorized_subset(
      tree = tree, X = X, rows = oob_idx,
      permute_var = permute_var, permuted_values = x_perm
    )
  }

  # The fitted tree structure remains fixed. The response distribution in each
  # leaf is based on either all original training observations (current estimator)
  # or, if explicitly requested, only the tree's in-bag observations.
  if (response_members == "inbag") {
    response_idx <- inbag
  } else {
    response_idx <- seq_len(n)
  }
  response_lids <- lids_base[response_idx]

  oob_by_leaf <- split(seq_along(oob_idx), oob_lids)
  response_by_leaf <- split(response_idx, response_lids)
  common_leaves <- intersect(names(oob_by_leaf), names(response_by_leaf))

  i_list <- vector("list", length(common_leaves))
  j_list <- vector("list", length(common_leaves))
  x_list <- vector("list", length(common_leaves))

  for (k in seq_along(common_leaves)) {
    lid <- common_leaves[k]
    row_pos <- oob_by_leaf[[lid]]
    target_idx <- oob_idx[row_pos]
    resp_idx <- response_by_leaf[[lid]]
    nL <- length(resp_idx)

    if (nL == 0L) next

    # For response_members = "all", exclude observation i from its own OOB
    # prediction whenever its original training covariates place it in the
    # terminal node being used for prediction. Under permutation, i may be
    # routed to a different leaf, in which case it is not present there anyway.
    if (response_members == "all") {
      denom_target <- nL - as.integer(target_idx %in% resp_idx)
    } else {
      # An OOB observation cannot be in this tree's in-bag sample.
      denom_target <- rep(nL, length(target_idx))
    }

    if (any(denom_target <= 0L)) {
      stop("Internal error: an OOB terminal node has no eligible response members.")
    }

    ii <- rep(row_pos, each = nL)
    jj <- rep(resp_idx, times = length(row_pos))
    target_rep <- rep(target_idx, each = nL)
    denom_rep <- rep(denom_target, each = nL)

    if (response_members == "all") {
      keep <- jj != target_rep
      ii <- ii[keep]
      jj <- jj[keep]
      denom_rep <- denom_rep[keep]
    }

    i_list[[k]] <- ii
    j_list[[k]] <- jj
    x_list[[k]] <- 1 / denom_rep
  }

  i_all <- unlist(i_list, use.names = FALSE)
  j_all <- unlist(j_list, use.names = FALSE)
  x_all <- unlist(x_list, use.names = FALSE)

  W <- Matrix::sparseMatrix(
    i = i_all, j = j_all, x = x_all,
    dims = c(length(oob_idx), n)
  )

  list(W = W, oob_idx = oob_idx)
}

# OOB permutation variable importance based on mean CRPS.
# For each variable, the fitted trees are kept fixed. The variable is permuted
# among OOB observations within each tree, and CRPS is calculated separately
# for each tree before averaging across trees. Under the current estimator, each
# OOB observation is predicted from all other training observations in its leaf.
# R controls the number of permutation repetitions.
variable_importance_crps.myRF <- function(object,
                                          R = 10, # number of repetitions per covariate
                                          variables = NULL, # which predictors receive permutation importance scores; if NULL, evaluate all covariates
                                          # which training outcomes may contribute to the terminal-node response distribution;
                                          # when omitted, the current estimator uses "all" with the target observation excluded
                                          response_members = c("inbag", "all"),
                                          seed = NULL, 
                                          sort = FALSE, # if sort=TRUE, output summary is ordered from largest to smallest importance
                                          verbose = FALSE # whether it prints progress messages
                                          ) {
  # Keep the public function signature unchanged. When omitted, use the
  # all-training-observations leaf distribution with leave-one-out exclusion.
  if (missing(response_members)) {
    response_members <- "all"
  } else {
    response_members <- match.arg(response_members)
  }
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || !is.finite(R) ||
      R < 1L || R != floor(R)) {
    stop("R must be a positive integer.")
  }
  R <- as.integer(R)

  xnames <- object$xnames
  p <- length(xnames)

  if (is.null(variables)) {
    var_idx <- seq_len(p)
  } else if (is.character(variables)) {
    var_idx <- match(variables, xnames)
    if (anyNA(var_idx)) stop("Some requested variables are not in object$xnames.")
  } else {
    var_idx <- as.integer(variables)
    if (anyNA(var_idx) || any(var_idx < 1L | var_idx > p)) {
      stop("Variable indices are out of range.")
    }
  }
  var_names <- xnames[var_idx]

  if (!is.null(seed)) set.seed(seed)

  ntree <- object$ntree
  X <- object$X
  y <- object$y

  # Cache the original training-data leaf memberships and the variables used
  # by each fitted tree. This avoids recomputing them for every permutation.
  all_lids <- vector("list", ntree)
  used_vars <- vector("list", ntree)
  for (t in seq_len(ntree)) {
    tree <- object$RF[[t]]
    all_lids[[t]] <- if (!is.null(tree$leaf_id_by_obs)) {
      tree$leaf_id_by_obs
    } else {
      .get_leaf_ids_vectorized(tree, X)
    }
    used_vars[[t]] <- .get_tree_used_vars(tree)
  }

  if (verbose) cat("Computing baseline tree-level OOB CRPS...\n")
  base_tree_mean <- rep(NA_real_, ntree)
  for (t in seq_len(ntree)) {
    tw <- .tree_oob_weights_for_vimp(
      object = object, tree_index = t,
      response_members = response_members,
      lids_base = all_lids[[t]],
      used_vars = used_vars[[t]]
    )
    if (length(tw$oob_idx) == 0L) next

    crps_t <- .calc_crps(
      as(tw$W, "RsparseMatrix"),
      y, y[tw$oob_idx]
    )
    base_tree_mean[t] <- mean(crps_t)
  }

  valid_trees <- is.finite(base_tree_mean)
  if (!any(valid_trees)) {
    stop("No trees have valid OOB observations for CRPS variable importance.")
  }
  base_mean <- mean(base_tree_mean[valid_trees])
  if (!is.finite(base_mean) || base_mean <= 0) {
    stop("Baseline mean OOB CRPS is not finite and positive.")
  }

  # Keep the per-observation OOB forest CRPS in the returned object for
  # compatibility with earlier versions. The variable-importance calculation
  # itself uses the tree-level losses above.
  base_crps <- predict_oob.myRF(object, what = "CRPS",
                                response_members = response_members)

  perm_mean <- matrix(NA_real_, nrow = R, ncol = length(var_idx),
                      dimnames = list(paste0("rep", seq_len(R)), var_names))

  for (r in seq_len(R)) {
    if (verbose) cat("Permutation repetition", r, "of", R, "\n")
    for (k in seq_along(var_idx)) {
      # We do NOT refit the forest!! The trees remain fixed.
      # We only permute one covariate among the OOB observations of each tree,
      # reroute those OOB observations through the already-grown tree, and
      # recalculate that tree's mean OOB CRPS.

      j <- var_idx[k]
      perm_tree_mean <- base_tree_mean

      for (t in which(valid_trees)) {
        # If the tree never split on this covariate, the permuted prediction
        # is identical to the baseline prediction for that tree.
        if (!(j %in% used_vars[[t]])) next

        tw <- .tree_oob_weights_for_vimp(
          object = object, tree_index = t,
          response_members = response_members,
          lids_base = all_lids[[t]],
          used_vars = used_vars[[t]],
          permute_var = j
        )

        crps_t <- .calc_crps(
          as(tw$W, "RsparseMatrix"),
          y, y[tw$oob_idx]
        )
        perm_tree_mean[t] <- mean(crps_t)
      }

      perm_mean[r, k] <- mean(perm_tree_mean[valid_trees])
      if (verbose) {
        # Because such variable importance takes a long time to compute,
        # we print progress messages for the user.
        cat("  ", var_names[k], ": mean tree-level OOB CRPS =",
            signif(perm_mean[r, k], 5), "\n")
      }
    }
  }

  # Average the permuted tree-level CRPS values over the R repetitions.
  perm_bar <- colMeans(perm_mean, na.rm = TRUE)
  vi_abs <- perm_bar - base_mean
  vi_pct <- 100 * vi_abs / base_mean # in percentage

  out <- data.frame(
    variable = var_names,
    baseline_mean_CRPS = rep(base_mean, length(var_idx)),
    permuted_mean_CRPS = as.numeric(perm_bar),
    VI_CRPS = as.numeric(vi_abs),
    VI_percent = as.numeric(vi_pct),
    row.names = NULL
  )

  if (sort) {
    out <- out[order(out$VI_percent, decreasing = TRUE), , drop = FALSE]
    rownames(out) <- NULL
  }

  structure(list(summary = out,
                 permuted_mean_CRPS = perm_mean,
                 baseline_CRPS = base_crps,
                 baseline_mean_CRPS = base_mean,
                 response_members = response_members,
                 R = R),
            class = "myRF_vimp_crps")
}

print.myRF_vimp_crps <- function(x, ...) {
  print(x$summary, row.names = FALSE, ...)
  invisible(x)
}


# Vectorized tree traversal.
# Instead of traversing one row at a time, send groups of row indices down the tree.
.get_leaf_ids_vectorized <- function(tree, X) {
  n <- nrow(X)
  leaf_ids <- integer(n)

  # Use an explicit stack rather than recursive superassignment. Each stack item
  # stores a node ID and the rows that have reached that node.
  stack <- list(list(node_id = 1L, idx = seq_len(n)))

  while (length(stack) > 0L) {
    current <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL

    node_id <- current$node_id
    idx <- current$idx
    if (length(idx) == 0L) next

    node <- tree$nodes[[node_id]]

    if (node$is_leaf) {
      # Save the terminal-node ID for all rows that reached this leaf.
      leaf_ids[idx] <- node_id
    } else {
      # Evaluate the stored split rule for this node.
      v <- node$covar
      xv <- X[[v]][idx]

      if (node$split_type == "numeric") {
        is_left <- xv <= node$split_point
      } else {
        # For categorical predictors, send rows left if their level is in the stored left set.
        is_left <- as.character(xv) %in% node$split_point
      }

      # Add right and left children to the stack. The order is not important for
      # the final leaf IDs.
      if (any(!is_left)) {
        stack[[length(stack) + 1L]] <- list(node_id = node$right, idx = idx[!is_left])
      }
      if (any(is_left)) {
        stack[[length(stack) + 1L]] <- list(node_id = node$left, idx = idx[is_left])
      }
    }
  }

  return(leaf_ids)
}

# Vectorized tree traversal for a subset of rows.
# Optionally override one predictor column with a permuted version. This is
# used for fixed-fitted-forest OOB permutation variable importance.
.get_leaf_ids_vectorized_subset <- function(tree, X, rows,
                                            permute_var = NULL,
                                            permuted_values = NULL) {
  if (length(rows) == 0L) return(integer(0))
  if (is.null(permute_var) != is.null(permuted_values)) {
    stop("Provide both 'permute_var' and 'permuted_values', or neither.")
  }

  n <- nrow(X)
  leaf_ids <- integer(n)
  rows <- as.integer(rows)

  # The stack stores the current tree node and the subset of requested rows that
  # reached that node. This avoids superassignment inside a nested recursive helper.
  stack <- list(list(node_id = 1L, idx = rows))

  while (length(stack) > 0L) {
    current <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL

    node_id <- current$node_id
    idx <- current$idx
    if (length(idx) == 0L) next

    node <- tree$nodes[[node_id]]

    if (node$is_leaf) {
      leaf_ids[idx] <- node_id
    } else {
      v <- node$covar
      if (!is.null(permute_var) & v == permute_var) {
        xv <- permuted_values[idx]
      } else {
        xv <- X[[v]][idx]
      }

      if (node$split_type == "numeric") {
        is_left <- xv <= node$split_point
      } else {
        is_left <- as.character(xv) %in% node$split_point
      }

      if (any(!is_left)) {
        stack[[length(stack) + 1L]] <- list(node_id = node$right, idx = idx[!is_left])
      }
      if (any(is_left)) {
        stack[[length(stack) + 1L]] <- list(node_id = node$left, idx = idx[is_left])
      }
    }
  }

  leaf_ids[rows]
}

# Calculates OOB weights using sparse matrices.
# Under the current estimator, an OOB observation is routed through a tree that
# did not use it for split selection, and its terminal-node distribution is formed
# from all other original training observations in that fitted leaf.
# response_members = "inbag" retains the historical strict in-bag alternative.
# If permute_var is not NULL, that covariate is permuted separately among the
# OOB observations of each tree before routing those OOB observations through
# the fixed fitted tree. The fitted tree and the stored training leaf memberships
# used to form terminal-node response distributions are not changed.
.oob_weights_fast <- function(object,
                              response_members = c("inbag", "all"),
                              permute_var = NULL) {
  response_members <- match.arg(response_members)
  n <- length(object$y)
  ntree <- object$ntree
  X <- object$X

  if (!is.null(permute_var)) {
    if (is.character(permute_var)) {
      permute_var <- match(permute_var, object$xnames)
      if (is.na(permute_var)) stop("'permute_var' is not in object$xnames.")
    }
    permute_var <- as.integer(permute_var)
    if (length(permute_var) != 1L || permute_var < 1L ||
        permute_var > length(object$xnames)) {
      stop("'permute_var' must be a single valid variable index or name.")
    }
  }

  # Leaf membership of the original training covariates under each fitted tree.
  # These memberships define the terminal-node response distributions.
  all_lids <- vector("list", ntree)
  used_vars <- vector("list", ntree)
  for (t in seq_len(ntree)) {
    tree <- object$RF[[t]]
    if (!is.null(tree$leaf_id_by_obs)) {
      all_lids[[t]] <- tree$leaf_id_by_obs
    } else {
      all_lids[[t]] <- .get_leaf_ids_vectorized(tree, X)
    }
    used_vars[[t]] <- .get_tree_used_vars(tree)
  }

  i_tree_list <- vector("list", ntree)
  j_tree_list <- vector("list", ntree)
  x_tree_list <- vector("list", ntree)
  denom <- integer(n)

  for (t in seq_len(ntree)) {
    tree <- object$RF[[t]]
    lids_base <- all_lids[[t]]
    inbag <- tree$inbag

    is_oob <- rep(TRUE, n)
    is_oob[unique(inbag)] <- FALSE
    oob_idx <- which(is_oob)
    if (length(oob_idx) == 0L) next

    if (is.null(permute_var) || !(permute_var %in% used_vars[[t]])) {
      # If the tree does not split on the permuted covariate, the OOB leaf
      # assignment is unchanged.
      oob_lids <- lids_base[oob_idx]
    } else {
      x_perm <- X[[permute_var]]
      x_perm[oob_idx] <- sample(x_perm[oob_idx], length(oob_idx), replace = FALSE)
      oob_lids <- .get_leaf_ids_vectorized_subset(
        tree = tree, X = X, rows = oob_idx,
        permute_var = permute_var, permuted_values = x_perm
      )
    }

    if (response_members == "inbag") {
      response_idx <- inbag
    } else {
      response_idx <- seq_len(n)
    }
    response_lids <- lids_base[response_idx]

    oob_by_leaf <- split(oob_idx, oob_lids)
    response_by_leaf <- split(response_idx, response_lids)
    common_leaves <- intersect(names(oob_by_leaf), names(response_by_leaf))

    i_leaf <- vector("list", length(common_leaves))
    j_leaf <- vector("list", length(common_leaves))
    x_leaf <- vector("list", length(common_leaves))

    for (k in seq_along(common_leaves)) {
      lid <- common_leaves[k]
      target_idx <- oob_by_leaf[[lid]]
      resp_idx <- response_by_leaf[[lid]]
      nL <- length(resp_idx)
      if (nL == 0L) next

      if (response_members == "all") {
        denom_target <- nL - as.integer(target_idx %in% resp_idx)
      } else {
        denom_target <- rep(nL, length(target_idx))
      }

      if (any(denom_target <= 0L)) {
        stop("Internal error: an OOB terminal node has no eligible response members.")
      }

      ii <- rep(target_idx, each = nL)
      jj <- rep(resp_idx, times = length(target_idx))
      target_rep <- rep(target_idx, each = nL)
      denom_rep <- rep(denom_target, each = nL)

      if (response_members == "all") {
        keep <- jj != target_rep
        ii <- ii[keep]
        jj <- jj[keep]
        denom_rep <- denom_rep[keep]
      }

      i_leaf[[k]] <- ii
      j_leaf[[k]] <- jj
      x_leaf[[k]] <- 1 / denom_rep
    }

    i_tree_list[[t]] <- unlist(i_leaf, use.names = FALSE)
    j_tree_list[[t]] <- unlist(j_leaf, use.names = FALSE)
    x_tree_list[[t]] <- unlist(x_leaf, use.names = FALSE)
    denom[oob_idx] <- denom[oob_idx] + 1L
  }

  i_all <- unlist(i_tree_list, use.names = FALSE)
  j_all <- unlist(j_tree_list, use.names = FALSE)
  x_all <- unlist(x_tree_list, use.names = FALSE)

  W_total <- Matrix::sparseMatrix(
    i = i_all, j = j_all, x = x_all,
    dims = c(n, n)
  )

  return(list(W = W_total, denom = denom))
}

# Add all-training leaf memberships to a tree object if not already present.
.ensure_leaf_membership_tree <- function(tree, X_all) {
  if (!is.null(tree$leaf_id_by_obs)) return(tree)
  
  tree$leaf_id_by_obs <- .get_leaf_ids_vectorized(tree, X_all)
  
  return(tree)
}


# Mean OOB pinball loss for a requested quantile level.
# This is a scalar diagnostic, not a splitting rule.
oob_pinball_loss.myRF <- function(object, tau = 0.5, na.rm = TRUE) {
  stopifnot(tau >= 0 & tau <= 1)
  y <- object$y
  Q <- predict_oob.myRF(object, what = "quantiles", probs = tau)
  q <- as.numeric(Q[, 1L])
  e <- y - q
  l <- .pinball_loss_vec(e, tau)
  if (na.rm) {
    l <- l[is.finite(l)]
    return(mean(l))
  } else {
    return(mean(l))
  }
}


# Extract the predictor indices used in splits in one tree.
# This allows permutation importance to skip re-routing trees that never used the permuted predictor.
.get_tree_used_vars <- function(tree) {
  vars <- integer()
  for (node in tree$nodes) {
    if (!node$is_leaf) {
      vars <- c(vars, node$covar)
    }
  }
  return(unique(vars))
}

