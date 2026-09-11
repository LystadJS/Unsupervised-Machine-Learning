# PCA / JSL Dark Academic v1.0.0 | 2026-09-10
# Dependencies: R >= 4.1; base and recommended packages only.
# Educational reference implementation around stats::prcomp(), not PCoA.
# All learned transformations are fitted on the supplied TRAINING rows only.

.pca_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
}

.pca_matrix <- function(x, features = NULL, min_rows = 1L) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("Input must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(trimws(colnames(x)))) || anyDuplicated(colnames(x))) {
    stop("Supply unique, nonempty feature names.", call. = FALSE)
  }
  if (is.null(features)) features <- colnames(x)
  if (!is.character(features) || !length(features) || anyNA(features) ||
      anyDuplicated(features) || any(!features %in% colnames(x))) {
    stop("features must be unique names present in the input.", call. = FALSE)
  }
  x <- x[, features, drop = FALSE]
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop("Selected features must be numeric; exclude IDs/outcomes/categories.",
         call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || is.complex(x) || nrow(x) < min_rows || !ncol(x)) {
    stop("Input must contain enough rows and real numeric features.", call. = FALSE)
  }
  if (any(is.infinite(x))) stop("Infinite values are not allowed.", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

.pca_k <- function(k, max_k, allow_zero = FALSE) {
  low <- if (allow_zero) 0L else 1L
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k != floor(k) || k < low || k > max_k) {
    stop("k must be an integer from ", low, " through ", max_k, ".", call. = FALSE)
  }
  as.integer(k)
}

.pca_object <- function(object) {
  if (!inherits(object, "jsl_pca")) stop("Expected a jsl_pca object.", call. = FALSE)
}

# Fit centered, ordinary PCA. scale_features=TRUE standardizes each feature.
# missing='median' learns one training median per feature; it is a simple policy,
# not a claim that median imputation is appropriate for every missingness process.
# constant='drop' removes only exactly zero-SD training features and records them.
# Near-constant features need substantive review and are NOT automatically dropped.
# k and variance_target are mutually exclusive. With neither, retain numerical rank.
# rank_tol is relative to the largest singular value; it is NOT a quality threshold.
pca_fit <- function(x, features = NULL, scale_features = TRUE,
                    missing = c("error", "median"),
                    constant = c("error", "drop"),
                    k = NULL, variance_target = NULL,
                    rank_tol = sqrt(.Machine$double.eps)) {
  .pca_flag(scale_features, "scale_features")
  missing <- match.arg(missing)
  constant <- match.arg(constant)
  if (!is.numeric(rank_tol) || length(rank_tol) != 1L ||
      !is.finite(rank_tol) || rank_tol <= 0 || rank_tol >= 1) {
    stop("rank_tol must be finite and strictly between zero and one.", call. = FALSE)
  }
  if (!is.null(k) && !is.null(variance_target)) {
    stop("Choose either k or variance_target, not both.", call. = FALSE)
  }
  x <- .pca_matrix(x, features, min_rows = 2L)
  feature_names <- colnames(x)
  if (any(colSums(!is.na(x)) == 0L)) {
    stop("A training feature is entirely missing; no median can be learned.", call. = FALSE)
  }
  na_count <- colSums(is.na(x))
  if (anyNA(x) && missing == "error") stop("Missing data: declare an explicit policy.", call. = FALSE)
  medians <- apply(x, 2L, stats::median, na.rm = TRUE)
  for (j in seq_len(ncol(x))) x[is.na(x[, j]), j] <- medians[j]
  sds <- apply(x, 2L, stats::sd)
  if (any(!is.finite(sds))) stop("Nonfinite standard deviation; rescale extreme values.", call. = FALSE)
  dropped <- names(sds)[sds == 0]
  if (length(dropped) && constant == "error") {
    stop("Constant training feature(s): ", paste(dropped, collapse = ", "), call. = FALSE)
  }
  keep <- setdiff(feature_names, dropped)
  if (!length(keep)) stop("No variable training features remain.", call. = FALSE)
  x <- x[, keep, drop = FALSE]
  fit <- stats::prcomp(x, center = TRUE, scale. = scale_features)
  eigenvalues <- fit$sdev^2
  total_variance <- sum(eigenvalues)
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop("PCA has no finite positive variance.", call. = FALSE)
  }
  rank <- sum(fit$sdev > rank_tol * max(fit$sdev))
  rank <- min(rank, nrow(x) - 1L, ncol(x))
  if (rank < 1L) stop("Numerical rank is zero.", call. = FALSE)
  proportion <- eigenvalues / total_variance
  cumulative <- pmin(1, cumsum(proportion))
  if (!is.null(variance_target)) {
    if (!is.numeric(variance_target) || length(variance_target) != 1L ||
        !is.finite(variance_target) || variance_target <= 0 || variance_target > 1) {
      stop("variance_target must be in (0, 1].", call. = FALSE)
    }
    candidates <- which(cumulative[seq_len(rank)] >= variance_target - 1e-12)
    if (!length(candidates)) {
      stop("Target exceeds variance retained at the chosen numerical rank tolerance.", call. = FALSE)
    }
    k <- candidates[1L]
  }
  if (is.null(k)) k <- rank
  k <- .pca_k(k, rank)
  object <- list(
    pca = fit, features = feature_names, kept_features = keep,
    dropped_features = dropped, medians = medians, missing = missing,
    scale_features = scale_features, training_sd = sds[keep],
    training_missing = na_count, n_train = nrow(x), rank = rank,
    rank_tol = rank_tol, k = k, variance_target = variance_target,
    spectrum = data.frame(component = seq_along(eigenvalues),
                          eigenvalue = eigenvalues, proportion = proportion,
                          cumulative = cumulative),
    version = "1.0.0")
  class(object) <- "jsl_pca"
  object
}

# Applies stored schema, training medians, and recorded feature removal.
# Extra columns are ignored. Missing required or duplicated column names stop.
# Row order is preserved; join external metadata by your own stable ID key.
pca_prepare <- function(object, newdata) {
  .pca_object(object)
  x <- .pca_matrix(newdata, object$features)
  if (anyNA(x) && object$missing == "error") {
    stop("New data contain missing values; fitted missing policy is 'error'.", call. = FALSE)
  }
  for (j in seq_len(ncol(x))) x[is.na(x[, j]), j] <- object$medians[j]
  x[, object$kept_features, drop = FALSE]
}

# Project new observations into the TRAINING frame; never refit on newdata.
pca_transform <- function(object, newdata, k = object$k) {
  .pca_object(object)
  k <- .pca_k(k, object$rank)
  x <- pca_prepare(object, newdata)
  stats::predict(object$pca, newdata = x)[, seq_len(k), drop = FALSE]
}

pca_scores <- function(object, k = object$k) {
  .pca_object(object)
  k <- .pca_k(k, object$rank)
  object$pca$x[, seq_len(k), drop = FALSE]
}

pca_loadings <- function(object, k = object$k) {
  .pca_object(object)
  k <- .pca_k(k, object$rank)
  object$pca$rotation[, seq_len(k), drop = FALSE]
}

# These are TRAINING variable-score correlations, not prcomp weight loadings.
# cor(X_j, T_l) = V_jl * sqrt(lambda_l) / sd(Z_j), with Z centered/scaled.
pca_variable_correlations <- function(object, k = object$k) {
  .pca_object(object)
  k <- .pca_k(k, object$rank)
  out <- sweep(pca_loadings(object, k), 2L, object$pca$sdev[seq_len(k)], "*")
  sd_z <- if (object$scale_features) rep(1, length(object$kept_features)) else object$training_sd
  sweep(out, 1L, sd_z, "/")
}

# Reconstruct kept columns only. k=0 gives the training centroid.
# original_units=FALSE returns the centered/scaled training feature space.
pca_reconstruct <- function(object, newdata = NULL, k = object$k,
                            original_units = TRUE) {
  .pca_object(object)
  .pca_flag(original_units, "original_units")
  k <- .pca_k(k, object$rank, allow_zero = TRUE)
  if (is.null(newdata)) {
    scores <- object$pca$x
  } else {
    x <- pca_prepare(object, newdata)
    scores <- stats::predict(object$pca, newdata = x)
  }
  if (k == 0L) {
    zhat <- matrix(0, nrow(scores), length(object$kept_features),
                   dimnames = list(rownames(scores), object$kept_features))
  } else {
    zhat <- scores[, seq_len(k), drop = FALSE] %*%
      t(object$pca$rotation[, seq_len(k), drop = FALSE])
  }
  if (original_units) {
    scale_vector <- if (object$scale_features) object$pca$scale else rep(1, ncol(zhat))
    zhat <- sweep(zhat, 2L, scale_vector, "*")
    zhat <- sweep(zhat, 2L, object$pca$center, "+")
  }
  zhat
}

# Error against preprocessed values (including any imputed cells); this is NOT
# an assessment of missing-value accuracy. Compare k only within the same frame.
pca_reconstruction_error <- function(object, data, k_grid = 0:object$rank) {
  .pca_object(object)
  if (!length(k_grid)) stop("k_grid must not be empty.", call. = FALSE)
  x <- pca_prepare(object, data)
  z <- scale(x, center = object$pca$center,
             scale = if (object$scale_features) object$pca$scale else FALSE)
  rows <- lapply(k_grid, function(k) {
    k <- .pca_k(k, object$rank, allow_zero = TRUE)
    zhat <- pca_reconstruct(object, data, k, original_units = FALSE)
    data.frame(k = k, mse = mean((z - zhat)^2),
               retained_training_variance = if (k == 0) 0 else object$spectrum$cumulative[k])
  })
  do.call(rbind, rows)
}

print.jsl_pca <- function(x, ...) {
  cat("Centered PCA |", x$n_train, "training rows |", length(x$kept_features), "features\n")
  cat("Scaled:", x$scale_features, "| Numerical rank:", x$rank, "| Retained:", x$k, "\n")
  cat("Retained training variance:", sprintf("%.2f%%", 100 * x$spectrum$cumulative[x$k]), "\n")
  if (length(x$dropped_features)) cat("Dropped constants:", paste(x$dropped_features, collapse = ", "), "\n")
  invisible(x)
}
