# Version 1.0.0 | 2026-09-10

validate_dissimilarity <- function(d, tol = 1e-10) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) ||
      tol <= 0 || tol >= 1) stop("tol must be one number strictly between 0 and 1.")
  if (!(inherits(d, "dist") || is.matrix(d)))
    stop("Supply a dist object or a full numeric square matrix, not a feature table.")
  D <- as.matrix(d)
  if (!is.numeric(D) || is.complex(D) || nrow(D) != ncol(D) || nrow(D) < 2L)
    stop("D must be a real numeric square matrix with at least two units.")
  if (any(!is.finite(D))) stop("D contains missing or non-finite distances.")
  if (any(D < 0)) stop("D contains negative distances; repair the input explicitly.")
  magnitude <- max(D)
  eps <- tol * magnitude
  if (max(abs(diag(D))) > eps) stop("The diagonal of D must be zero.")
  if (max(abs(D - t(D))) > eps) stop("D is asymmetric; do not discard one triangle.")
  rn <- rownames(D); cn <- colnames(D)
  if (xor(is.null(rn), is.null(cn))) stop("Supply both row and column IDs, or neither.")
  if (!is.null(rn)) {
    if (!identical(rn, cn)) stop("Row and column IDs must match in the same order.")
    if (anyNA(rn) || any(!nzchar(trimws(rn))) || anyDuplicated(rn))
      stop("Unit IDs must be unique, nonmissing, and nonempty.")
  } else {
    rn <- sprintf("unit_%03d", seq_len(nrow(D)))
    dimnames(D) <- list(rn, rn)
  }
  adjustment <- 
    max(abs(D - t(D)), abs(diag(D)))
  D <- (D / 2 + t(D) / 2)
  diag(D) <- 0
  attr(D, "validation_adjustment_max") <- adjustment
  D
}

pcoa_gram <- function(D) {
  Q <- D^2
  if (any(!is.finite(Q))) stop("Squared distances overflow; rescale distance units.")
  B <- -0.5 * (sweep(sweep(Q, 1L, rowMeans(Q), "-"),
                     2L, colMeans(Q), "-") + mean(Q))
  (B / 2 + t(B) / 2)
}

pcoa_fit <- function(d, k = 2L, correction = c("none", "lingoes", "cailliez"),
                     tol = 1e-10) {
  correction <- match.arg(correction)
  D <- validate_dissimilarity(d, tol)
  n <- nrow(D)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != as.integer(k) ||
      k < 1L || k >= n) stop("k must be an integer from 1 through n - 1.")
  unit_scale <- max(D)
  if (unit_scale == 0) stop("All distances are zero: no positive axes or inertia shares exist.")
  if (!is.finite(unit_scale^2) || unit_scale^2 == 0)
    stop("Distance units are outside the supported numeric range; rescale D.")

  D0 <- D / unit_scale
  original <- eigen(pcoa_gram(D0), symmetric = TRUE)
  cutoff0 <- tol * max(abs(original$values))
  neg0 <- original$values < -cutoff0
  working <- D0
  c_normalized <- 0
  if (any(neg0) && correction == "lingoes") {
    c_normalized <- -min(original$values)
    working <- sqrt(D0^2 + 2 * c_normalized)
    diag(working) <- 0
  }
  if (any(neg0) && correction == "cailliez") {
    corrected <- stats::cmdscale(stats::as.dist(D0), k = 1L,
                                eig = TRUE, add = TRUE)
    c_normalized <- max(0, corrected$ac)
    working <- D0 + c_normalized
    diag(working) <- 0
  }
  spectral <- eigen(pcoa_gram(working), symmetric = TRUE)
  cutoff <- tol * max(abs(spectral$values))
  positive <- which(spectral$values > cutoff)
  negative <- which(spectral$values < -cutoff)
  if (!length(positive)) stop("No positive coordinate axes survive the numerical tolerance.")
  if (correction != "none" && length(negative))
    stop("Correction left substantive negative eigenvalues; inspect scale and tolerance.")
  if (any(neg0) && correction == "none")
    warning("Non-Euclidean input: positive-axis coordinates omit negative components.",
            call. = FALSE)
  used_k <- min(as.integer(k), length(positive))
  if (used_k < k) warning("Fewer positive axes than requested; returning ", used_k, ".",
                         call. = FALSE)
  idx <- positive[seq_len(used_k)]
  Z <- sweep(spectral$vectors[, idx, drop = FALSE], 2L,
             sqrt(spectral$values[idx]), "*") * unit_scale

  for (j in seq_len(ncol(Z))) {
    anchor <- which.max(abs(Z[, j]))
    if (Z[anchor, j] < 0) Z[, j] <- -Z[, j]
  }
  dimnames(Z) <- list(rownames(D), paste0("PCoA", seq_len(used_k)))
  eig_work <- spectral$values * unit_scale^2
  eig_orig <- original$values * unit_scale^2
  positive_sum <- sum(spectral$values[positive])
  abs_sum <- sum(abs(spectral$values[abs(spectral$values) > cutoff]))
  orig_abs_sum <- sum(abs(original$values[abs(original$values) > cutoff0]))
  working_distance <- working * unit_scale
  dimnames(working_distance) <- dimnames(D)
  additive_constant <- c_normalized * if (correction == "lingoes") unit_scale^2 else unit_scale
  shares <- ifelse(spectral$values > cutoff, spectral$values / positive_sum, 0)
  result <- list(
    points = Z,
    eigenvalues = data.frame(axis = seq_len(n), original = eig_orig,
      working = eig_work, positive_share = shares,
      cumulative_positive_share = cumsum(shares)),
    original_distance = D, working_distance = working_distance,
    correction_requested = correction,
    correction_applied = if (c_normalized > 0) correction else "none",
    additive_constant = additive_constant,
    constant_units = if (correction == "lingoes") "squared distance units (sqrt(d^2 + 2c))" else "distance units (d + c)",
    k_requested = as.integer(k), k_returned = used_k,
    positive_rank = length(positive), negative_count_original = sum(neg0),
    negative_count_working = length(negative),
    negative_fraction_original = sum(abs(original$values[neg0])) / orig_abs_sum,
    gof_absolute = sum(spectral$values[idx]) / abs_sum,
    gof_positive = sum(spectral$values[idx]) / positive_sum,
    axis_share_basis = if (c_normalized > 0) "corrected positive inertia" else "original positive inertia",
    tolerance_relative = tol,
    eigenvalue_cutoff_original = cutoff0 * unit_scale^2,
    eigenvalue_cutoff_working = cutoff * unit_scale^2,
    validation_adjustment_max = attr(D, "validation_adjustment_max")
  )
  class(result) <- "jsl_pcoa"
  result
}

pcoa_diagnostics <- function(fit) {
  if (!inherits(fit, "jsl_pcoa")) stop("fit must come from pcoa_fit().")
  keep <- lower.tri(fit$original_distance)
  dhat <- as.matrix(stats::dist(fit$points))[keep]
  compare <- function(D, target) {
    d <- D[keep]
    scale <- max(d)
    residual <- (dhat - d) / scale
    rho <- if (stats::sd(d) > 0 && stats::sd(dhat) > 0)
      stats::cor(d, dhat, method = "spearman") else NA_real_
    data.frame(target = target, pairs = length(d),
      normalized_distance_residual = sqrt(sum(residual^2) / sum((d / scale)^2)),
      max_absolute_error = max(abs(dhat - d)), spearman_rho = rho)
  }
  rbind(compare(fit$original_distance, "original"),
        compare(fit$working_distance, "working"))
}

pcoa_scores <- function(fit, metadata = NULL, id_col = "id") {
  if (!inherits(fit, "jsl_pcoa")) stop("fit must come from pcoa_fit().")
  scores <- data.frame(id = rownames(fit$points), fit$points, row.names = NULL,
                       check.names = FALSE)
  if (is.null(metadata)) return(scores)
  if (!is.data.frame(metadata) || !id_col %in% names(metadata))
    stop("metadata must be a data frame containing id_col.")
  ids <- as.character(metadata[[id_col]])
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
    stop("Metadata IDs must be unique, nonmissing, and nonempty.")
  if (!setequal(ids, scores$id)) stop("Metadata IDs must exactly match ordination IDs.")
  other <- setdiff(names(metadata), id_col)
  if (any(other %in% names(scores))) stop("Metadata names collide with coordinate names.")
  cbind(scores, metadata[match(scores$id, ids), other, drop = FALSE], row.names = NULL)
}

pcoa_axis_labels <- function(fit) {
  if (!inherits(fit, "jsl_pcoa")) 
    stop("fit must come from pcoa_fit().")
  s <- 
    head(fit$eigenvalues$positive_share, fit$k_returned)
  sprintf("PCoA %d (%.1f%% %s)", seq_along(s), 100 * s, fit$axis_share_basis)
}

# PLOTTING:
plot_pcoa <- function(fit, metadata = NULL, id_col = "id", group_col = NULL, dark = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_pcoa() requires ggplot2. Core fitting and diagnostics do not.")
  if (!is.logical(dark) || length(dark) != 1L || is.na(dark)) stop("dark must be TRUE or FALSE.")
  dat <- 
    pcoa_scores(fit, metadata, id_col)
  if (fit$k_returned < 2L) stop("A two-axis plot requires two returned positive axes.")
  if (!is.null(group_col) && !group_col %in% names(dat)) stop("Unknown group_col.")
  labels <- 
    pcoa_axis_labels(fit)
  g <- 
    ggplot2::ggplot(dat, ggplot2::aes(x = PCoA1, y = PCoA2))
  if (is.null(group_col)) {
    g <- 
      g + ggplot2::geom_point(size = 2.6, colour = if (dark) "#EEE5D9" else "#312632")
  } else {
    dat$.group <- 
        factor(dat[[group_col]])
    if (anyNA(dat$.group)) stop("Group labels contain missing values.")
    if (nlevels(dat$.group) > 4L) stop("The compact JSL palette supports up to four groups.")
    colours <- if (dark) c("#C8BCC7", "#A6425D", "#9B80AE", "#92A9BA") else
      c("#312632", "#8D1732", "#65457E", "#4D6B7B")
    g <- 
        ggplot(
          dat, 
          aes(
            x     = PCoA1, 
            y     = PCoA2, 
            color = .group, 
            shape = .group
          )
        ) +
      geom_point(
        size = 2.5
      ) +
      scale_colour_manual(
        values = colors
      ) +
      labs(
        colour = group_col, 
        shape  = group_col
      )
  }
  bg <- if (dark) "#09080B" else "#FFFFFF"
  fg <- if (dark) "#C8BCC7" else "#312632"
  g + coord_fixed(ratio = 1) +
    labs(
      x        = labels[1], 
      y        = labels[2],
      title    = "Proximity summarizes the chosen dissimilarity",
      subtitle = paste("Correction:", fit$correction_applied),
      caption  = "A descriptive map; group separation is not a significance test."
    ) +
    theme_minimal(
      base_size   = 11, 
      base_family = "sans"
    ) +
    theme(
      panel.grid        = element_blank(),
      plot.background   = element_rect(fill = bg, color = NA),
      panel.background  = element_rect(fill = bg, color = NA),
      text              = element_text(color = fg),
      axis.text         = element_text(color = fg),
      legend.background = element_rect(fill = bg, color = NA),
      legend.key        = element_rect(fill = bg, color = NA))
}
