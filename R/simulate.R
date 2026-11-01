# Simulation engine for the matrix-normal LRCP/LRCQ model.
#
# Scheme (slides "Simulation scheme: matrix normal"):
#   B (m x q) ~ MN(0, U, V),  U = W^{1/2} U0 W^{1/2},  V = Hg / M
#   so that Cov(b_.a, b_.b) = (h_ab / M) U and Sigma_beta = avg(h2)/M * U.
#   Z = R B diag(sqrt(n)) + E,   E ~ MN(0, R, C)  (C = intercept matrix).

#' Symmetric square root factor of a PSD matrix
#'
#' Returns L with \code{L \%*\% t(L) == S}, using Cholesky when possible and an
#' eigen decomposition (negative eigenvalues floored at zero) otherwise.
#' @param S Symmetric PSD matrix.
#' @keywords internal
psd_factor <- function(S) {
  S <- as.matrix(S)
  L <- tryCatch(t(chol(S)), error = function(e) NULL)
  if (!is.null(L)) return(L)
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
}

#' Simulate genetic-effect enrichment
#'
#' Draws \code{w} as Bernoulli(\code{prop_nonzero}) times a positive draw and
#' rescales to mean 1, so \code{sum(w) == m} as required by the definition
#' \eqn{w = diag(\Sigma_\beta)/avg(diag(\Sigma_\beta))}.
#'
#' @param m Number of SNPs.
#' @param prop_nonzero Proportion of SNPs with non-zero enrichment.
#' @param dist Distribution of non-zero values: log-normal, exponential, or a
#'   point mass (all non-zero SNPs equal).
#' @param sdlog Log-scale SD for \code{dist = "lnorm"}.
#' @param idx Optional indices of non-zero SNPs (overrides random selection).
#' @return Length-m vector with mean 1.
#' @export
make_w <- function(m, prop_nonzero = 0.1, dist = c("lnorm", "exp", "point"),
                   sdlog = 1, idx = NULL) {
  dist <- match.arg(dist)
  if (is.null(idx)) {
    k <- max(1L, round(prop_nonzero * m))
    idx <- sort(sample.int(m, k))
  }
  w <- numeric(m)
  w[idx] <- switch(dist,
    lnorm = stats::rlnorm(length(idx), 0, sdlog),
    exp = stats::rexp(length(idx)),
    point = rep(1, length(idx)))
  w * m / sum(w)
}

#' Simulate a genetic-effect correlation matrix
#'
#' @param m Number of SNPs.
#' @param type \code{"identity"}: no correlation. \code{"cluster"}: SNPs in
#'   \code{clusters} share an exchangeable correlation \code{rho} within each
#'   cluster, with optional random sign flips (still positive definite).
#'   \code{"pairs"}: explicit pairs with given correlations.
#' @param clusters List of integer vectors (SNP indices) for \code{"cluster"}.
#' @param rho Within-cluster correlation (scalar or one per cluster).
#' @param signs Randomly flip the sign of SNPs within clusters.
#' @param pairs Two-column matrix of SNP index pairs for \code{"pairs"}.
#' @param pair_rho Correlation for each pair.
#' @return m-by-m correlation matrix (class \code{Matrix} when sparse).
#' @export
make_Rb <- function(m, type = c("identity", "cluster", "pairs"),
                    clusters = NULL, rho = 0.5, signs = FALSE,
                    pairs = NULL, pair_rho = 0.5) {
  type <- match.arg(type)
  Rb <- Matrix::Diagonal(m)
  Rb <- methods::as(Rb, "generalMatrix")
  if (type == "cluster") {
    rho <- rep_len(rho, length(clusters))
    for (c in seq_along(clusters)) {
      ix <- clusters[[c]]
      k <- length(ix)
      if (k < 2) next
      if (rho[c] <= -1 / (k - 1) || rho[c] >= 1)
        stop("cluster rho outside the positive-definite range")
      s <- if (signs) sample(c(-1, 1), k, replace = TRUE) else rep(1, k)
      blk <- rho[c] * outer(s, s)
      diag(blk) <- 1
      Rb[ix, ix] <- blk
    }
  } else if (type == "pairs") {
    pairs <- matrix(pairs, ncol = 2)
    pair_rho <- rep_len(pair_rho, nrow(pairs))
    for (p in seq_len(nrow(pairs))) {
      Rb[pairs[p, 1], pairs[p, 2]] <- pair_rho[p]
      Rb[pairs[p, 2], pairs[p, 1]] <- pair_rho[p]
    }
    ev <- min(eigen(as.matrix(Rb), symmetric = TRUE, only.values = TRUE)$values)
    if (ev <= 0) stop("pairs give a non positive-definite Rb; ",
                      "keep each SNP in at most one pair or lower pair_rho")
  }
  Matrix::forceSymmetric(Rb)
}

#' Simulate trait-level genetic covariance
#'
#' @param q Number of traits.
#' @param h2 Heritabilities (length q); if \code{NULL}, drawn from
#'   Beta(\code{shape1}, \code{shape2}) scaled into \code{h2_range}.
#' @param type Trait genetic-correlation structure: none, exchangeable
#'   (\code{"cs"}), or block clusters of traits.
#' @param rg Genetic correlation within the structure.
#' @param n_clusters Number of trait clusters for \code{"cluster"}.
#' @param shape1,shape2 Beta distribution shapes for drawing heritabilities.
#' @param h2_range Range the Beta draws are scaled into.
#' @return List with \code{h2}, \code{Rg} and \code{gcov}
#'   (\eqn{h_{ab} = r_{g,ab} h_a h_b}).
#' @export
make_gcov <- function(q, h2 = NULL, type = c("identity", "cs", "cluster"),
                      rg = 0.3, n_clusters = 3, shape1 = 2, shape2 = 5,
                      h2_range = c(0.02, 0.6)) {
  type <- match.arg(type)
  if (is.null(h2)) {
    h2 <- h2_range[1] + diff(h2_range) * stats::rbeta(q, shape1, shape2)
  }
  Rg <- diag(q)
  if (type == "cs") {
    Rg[] <- rg
    diag(Rg) <- 1
  } else if (type == "cluster") {
    lab <- rep_len(seq_len(n_clusters), q)
    Rg <- outer(lab, lab, "==") * rg
    diag(Rg) <- 1
  }
  sh <- sqrt(h2)
  list(h2 = h2, Rg = Rg, gcov = Rg * outer(sh, sh))
}

#' Inter-GWAS intercept matrix
#'
#' \eqn{C_{ab} = r_{ab} o_{ab}} where \eqn{r_{ab}} is the phenotypic
#' (residual) correlation and \eqn{o_{ab} = N_s/\sqrt{n_a n_b}} the sample
#' overlap; the diagonal is the univariate intercept (1 without
#' stratification).
#'
#' @param q Number of traits.
#' @param overlap Scalar or q-by-q overlap fractions.
#' @param rp Scalar or q-by-q phenotypic correlations.
#' @param diag_intercept Univariate intercepts (scalar or length q).
#' @export
make_intercept <- function(q, overlap = 0, rp = 0, diag_intercept = 1) {
  O <- if (length(overlap) == 1) matrix(overlap, q, q) else overlap
  P <- if (length(rp) == 1) matrix(rp, q, q) else rp
  C <- O * P
  diag(C) <- rep_len(diag_intercept, q)
  C
}

#' Simulate the genetic-effect matrix B
#'
#' Draws \eqn{B \sim MN(0, U, V)} with \eqn{U = W^{1/2} R_\beta W^{1/2}} and
#' \eqn{V = H_g / M}.
#'
#' @param w Enrichment (length m).
#' @param Rb Genetic-effect correlation (m by m) or \code{NULL}.
#' @param gcov q-by-q genetic covariance (or vector of heritabilities).
#' @param M Number of SNPs heritability is spread over (default m).
#' @return m-by-q matrix.
#' @export
simulate_B <- function(w, Rb = NULL, gcov, M = length(w)) {
  m <- length(w)
  if (is.null(dim(gcov))) gcov <- diag(gcov, length(gcov))
  q <- nrow(gcov)
  X <- matrix(stats::rnorm(m * q), m, q)
  if (!is.null(Rb)) {
    nz <- which(w > 0)
    Lb <- psd_factor(as.matrix(Rb)[nz, nz, drop = FALSE])
    Xn <- matrix(0, m, q)
    Xn[nz, ] <- Lb %*% X[nz, , drop = FALSE]
    X <- Xn
  }
  Lv <- psd_factor(gcov / M)
  sqrt(w) * X %*% t(Lv)
}

#' Apply an LD matrix given as a matrix or a list of diagonal blocks
#' @keywords internal
ld_apply <- function(ld, X, fun = function(R, X) R %*% X) {
  if (is.list(ld) && !inherits(ld, "Matrix")) {
    sizes <- vapply(ld, nrow, 1L)
    end <- cumsum(sizes)
    start <- end - sizes + 1
    out <- vector("list", length(ld))
    for (b in seq_along(ld)) {
      out[[b]] <- fun(ld[[b]], X[start[b]:end[b], , drop = FALSE])
    }
    do.call(rbind, out)
  } else {
    fun(ld, X)
  }
}

#' Simulate GWAS Z-scores from genetic effects
#'
#' \eqn{Z = R B\, diag(\sqrt n) + E}, \eqn{E \sim MN(0, R, C)}. This is the
#' summary-level generative model of chapter 6 (\eqn{\hat A = RB + U}).
#'
#' @param ld LD correlation matrix, or a list of LD blocks (block diagonal).
#' @param B m-by-q genetic effects.
#' @param n Length-q sample sizes.
#' @param intercept q-by-q intercept matrix (see [make_intercept()]).
#' @return m-by-q Z-score matrix.
#' @export
simulate_Z <- function(ld, B, n, intercept = diag(length(n))) {
  m <- nrow(B)
  q <- ncol(B)
  signal <- ld_apply(ld, B) %*% diag(sqrt(n), q)
  G <- matrix(stats::rnorm(m * q), m, q)
  noise <- ld_apply(ld, G, function(R, X) psd_factor(R) %*% X)
  noise <- noise %*% t(psd_factor(intercept))
  as.matrix(signal + noise)
}

#' Simulate a full LRCP/LRCQ data set
#'
#' Convenience wrapper: draws B and Z and returns them with the true
#' parameters.
#'
#' @param ld LD matrix or list of blocks.
#' @param n Sample sizes (length q, or scalar recycled to \code{q}).
#' @param q Number of traits (used when \code{gcov} is \code{NULL}).
#' @param w Enrichment; drawn with [make_w()] when \code{NULL}.
#' @param Rb Genetic-effect correlation (\code{NULL} for none).
#' @param gcov Genetic covariance; drawn with [make_gcov()] when \code{NULL}.
#' @param intercept Intercept matrix; identity when \code{NULL}.
#' @param M Number of SNPs heritability is spread over (default m).
#' @param ... Passed to [make_w()].
#' @return List with \code{Z}, \code{B} and the true parameters.
#' @export
simulate_lrcpq <- function(ld, n, q = length(n), w = NULL, Rb = NULL,
                           gcov = NULL, intercept = NULL, M = NULL, ...) {
  m <- if (is.list(ld) && !inherits(ld, "Matrix")) sum(vapply(ld, nrow, 1L)) else nrow(ld)
  n <- rep_len(n, q)
  if (is.null(w)) w <- make_w(m, ...)
  if (is.null(gcov)) gcov <- make_gcov(q)$gcov
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(M)) M <- m
  B <- simulate_B(w, Rb, gcov, M)
  Z <- simulate_Z(ld, B, n, intercept)
  list(Z = Z, B = B, w = w, Rb = Rb, gcov = gcov, h2 = diag(gcov), n = n,
       intercept = intercept, M = M)
}

#' Toy LD: block-diagonal AR(1) correlation
#'
#' @param block_sizes Integer vector of block sizes.
#' @param rho AR(1) correlation between adjacent SNPs within a block.
#' @return List of LD blocks.
#' @export
make_ld_ar1 <- function(block_sizes, rho = 0.5) {
  lapply(block_sizes, function(k) rho^abs(outer(seq_len(k), seq_len(k), "-")))
}

#' Bind a list of LD blocks into one (sparse) block-diagonal matrix
#' @param blocks List of LD blocks.
#' @export
ld_bdiag <- function(blocks) Matrix::bdiag(blocks)
