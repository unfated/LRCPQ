#' Genetic-effect covariance kernel
#'
#' Builds the m-by-m kernel \eqn{S_w = W^{1/2} R_\beta W^{1/2}} so that the
#' SNP-wise genetic-effect covariance for traits a and b is
#' \eqn{\Sigma_{ab} = (h_{ab}/M) S_w}. Its diagonal is the enrichment vector
#' \code{w}; its off-diagonal holds \eqn{\rho_{ij}\sqrt{w_i w_j}}.
#'
#' @param w Length-m enrichment vector (mean 1 genome-wide).
#' @param Rb m-by-m genetic-effect correlation matrix, or \code{NULL} for the
#'   identity (no genetic-effect correlation).
#' @param pair_sum How off-diagonal pairs enter the moment equations.
#'   \code{"once"} (default) counts each unordered SNP pair once, which is what
#'   expanding \eqn{r_k' \Sigma r_l} gives. \code{"chapter"} reproduces the
#'   thesis chapter 6 equations, which sum the cross term over all ordered pairs
#'   and so count each pair twice (off-diagonal doubled).
#' @return An m-by-m matrix (dense or \code{Matrix}).
#' @export
sigma_w <- function(w, Rb = NULL, pair_sum = c("once", "chapter")) {
  pair_sum <- match.arg(pair_sum)
  m <- length(w)
  if (any(w < 0)) stop("w must be non-negative to build a covariance kernel")
  if (is.null(Rb)) return(diag(w, m))
  if (!all(dim(Rb) == m)) stop("Rb must be m-by-m with m = length(w)")
  sw <- sqrt(w)
  S <- as.matrix(Rb) * outer(sw, sw)
  if (pair_sum == "chapter") {
    d <- diag(S)
    S <- 2 * S
    diag(S) <- d
  }
  diag(S) <- w
  S
}

#' Per-trait scaling constants
#'
#' \eqn{c_{ab} = \sqrt{n_a n_b}\, h_{ab} / M}, the coefficient multiplying the LD
#' design in the moment equations. For a = b this is \eqn{s_a = n_a h^2_a / M}
#' of the fast OLS derivation.
#'
#' @param n Length-q GWAS sample sizes.
#' @param gcov q-by-q genetic covariance matrix (heritabilities on the
#'   diagonal, co-heritabilities off it), or a length-q vector of
#'   heritabilities.
#' @param M Number of SNPs over which heritability is defined (genome-wide).
#' @return q-by-q matrix of \eqn{c_{ab}}.
#' @export
scale_matrix <- function(n, gcov, M) {
  if (is.null(dim(gcov))) gcov <- diag(gcov, length(gcov))
  sn <- sqrt(n)
  outer(sn, sn) * gcov / M
}

#' Expected Z-score cross-product matrix
#'
#' Matrix form of the LRCP/LRCQ moment equation (chapter 6, eq. 6.2.4.3):
#' \deqn{E[z_{\cdot a} z_{\cdot b}'] = c_{ab} R S_w R + C_{ab} R,}
#' with \eqn{c_{ab} = \sqrt{n_a n_b} h_{ab}/M} and \eqn{C_{ab}} the
#' inter-GWAS intercept (\eqn{r_{ab} o_{ab}} off the diagonal, the univariate
#' LDSC intercept on it).
#'
#' @param R m-by-m LD correlation matrix.
#' @param w,Rb,pair_sum See [sigma_w()].
#' @param n,gcov,M See [scale_matrix()].
#' @param intercept q-by-q intercept matrix (default identity).
#' @param a,b Trait indices.
#' @return m-by-m matrix.
#' @export
expected_zz <- function(R, w, Rb = NULL, n, gcov, M, intercept = NULL,
                        a = 1, b = a, pair_sum = c("once", "chapter")) {
  q <- length(n)
  if (is.null(intercept)) intercept <- diag(q)
  cab <- scale_matrix(n, gcov, M)[a, b]
  R <- as.matrix(R)
  cab * (R %*% sigma_w(w, Rb, pair_sum) %*% R) + intercept[a, b] * R
}

#' Expected squared Z-scores (LRCQ moment)
#'
#' \eqn{E[z_{ka}^2] = s_a \{(R S_w R)_{kk}\} + C_{aa}}. With \code{Rb = NULL}
#' this reduces to \eqn{s_a \sum_s w_s r_{ks}^2 + C_{aa}}, the chi-square-only
#' LRCQ model.
#'
#' @inheritParams expected_zz
#' @return m-by-q matrix.
#' @export
expected_z2 <- function(R, w, Rb = NULL, n, gcov, M, intercept = NULL,
                        pair_sum = c("once", "chapter")) {
  q <- length(n)
  if (is.null(intercept)) intercept <- diag(q)
  h2 <- if (is.null(dim(gcov))) gcov else diag(gcov)
  s <- n * h2 / M
  R <- as.matrix(R)
  if (is.null(Rb)) {
    g <- as.vector((R^2) %*% w)
  } else {
    g <- rowSums((R %*% sigma_w(w, Rb, pair_sum)) * R)
  }
  outer(g, s) + matrix(diag(intercept), length(g), q, byrow = TRUE)
}

#' Literal element-wise chapter 6 moment equation
#'
#' Evaluates \eqn{E[z_{ka} z_{lb}]} term by term as written in chapter 6
#' (\eqn{i \ne j} in the double sum). Used to test the matrix form; slow.
#'
#' @param k,l SNP indices.
#' @param w Enrichment vector.
#' @param Rb Genetic-effect correlation matrix.
#' @inheritParams expected_zz
#' @param pair_factor Multiplier on the double sum: 1 reproduces the chapter
#'   equation literally (each pair counted twice), 0.5 counts each pair once.
#' @export
expected_cross_element <- function(k, l, a, b, R, w, Rb, n, gcov, M,
                                   intercept = NULL, pair_factor = 1) {
  q <- length(n)
  if (is.null(intercept)) intercept <- diag(q)
  cab <- scale_matrix(n, gcov, M)[a, b]
  R <- as.matrix(R)
  m <- length(w)
  first <- sum(w * R[k, ] * R[l, ])
  second <- 0
  for (i in seq_len(m)) for (j in seq_len(m)) {
    if (i == j) next
    second <- second + Rb[i, j] * sqrt(w[i] * w[j]) *
      (R[i, l] * R[j, k] + R[i, k] * R[j, l])
  }
  cab * (first + pair_factor * second) + intercept[a, b] * R[k, l]
}

#' Expected cross-trait covariance of one SNP's Z-scores
#'
#' \eqn{M_k = c\,(R S_w R)_{kk} + C}, the q-by-q covariance of
#' \eqn{(z_{k1}, \dots, z_{kq})}. By Isserlis' theorem
#' \eqn{Cov(z_{ka}^2, z_{kb}^2) = 2 M_{k,ab}^2} and, for SNPs without LD,
#' \eqn{Cov(z_{ka} z_{la}, z_{kb} z_{lb}) = M_{k,ab} M_{l,ab}}. The effective
#' number of traits for a distal pair is
#' \eqn{q_{eff} = tr(M_k) tr(M_l) / tr(M_k M_l)}.
#'
#' @param k SNP index (or vector of indices) within \code{R}.
#' @inheritParams expected_zz
#' @return A q-by-q matrix, or a list of them when \code{k} has length > 1.
#' @export
expected_trait_cov <- function(k, R, w, Rb = NULL, n, gcov, M, intercept = NULL,
                               pair_sum = c("once", "chapter")) {
  q <- length(n)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  cmat <- scale_matrix(n, gcov, M)
  R <- as.matrix(R)
  g <- rowSums((R[k, , drop = FALSE] %*% sigma_w(w, Rb, pair_sum)) * R[k, , drop = FALSE])
  out <- lapply(g, function(gk) gk * cmat + intercept)
  if (length(k) == 1) out[[1]] else out
}

#' Effective number of traits for a SNP pair
#'
#' A single-number summary that is exact for null standard errors only; for
#' power use [effective_traits()] (\eqn{q_P}, \eqn{q_G}, \eqn{q_{PG}}).
#' @param Mk,Ml Trait covariance matrices from [expected_trait_cov()].
#' @export
q_eff <- function(Mk, Ml = Mk) sum(diag(Mk)) * sum(diag(Ml)) / sum(Mk * Ml)
