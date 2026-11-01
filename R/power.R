# Power and sample size (theory supplement S9).

#' Effective numbers of traits
#'
#' \eqn{q_P = q^2/\sum\Gamma_{ab}^2} (sample overlap times phenotypic
#' correlation), \eqn{q_G = q^2/\sum r_{g,ab}^2} and
#' \eqn{q_{PG} = q^2/\sum \Gamma_{ab} r_{g,ab}}. Null standard errors depend
#' only on \eqn{q_P}; power at large enrichment is governed by \eqn{q_G}
#' (theory supplement S9.1). A single \eqn{q_{eff}} is only a null summary.
#'
#' @param intercept q-by-q intercept matrix \eqn{\Gamma} (scaled to unit
#'   diagonal internally).
#' @param rg q-by-q genetic correlation matrix (or a genetic covariance,
#'   converted to correlation).
#' @export
effective_traits <- function(intercept, rg) {
  G <- stats::cov2cor(as.matrix(intercept))
  r <- stats::cov2cor(as.matrix(rg))
  q <- nrow(G)
  c(q_P = q^2 / sum(G^2), q_G = q^2 / sum(r^2), q_PG = q^2 / sum(G * r))
}

#' Standard error of a clump enrichment estimate (eq. 9.1)
#'
#' @param w Enrichment (clump total).
#' @param n,h2 Typical GWAS sample size and heritability.
#' @param M Number of SNPs fixing the per-SNP scale (\eqn{s = n h^2/M}).
#' @param qe Effective trait numbers from [effective_traits()] (named
#'   vector), or one number used for all three.
#' @param kappa \eqn{(D^{-1})_{ii}} for tags (1.4 after \eqn{r^2 < 0.5}
#'   pruning).
#' @export
se_lrcq <- function(w, n, h2, M = 1.2e6, qe, kappa = 1.4) {
  qe <- qe_vec(qe)
  s <- n * h2 / M
  sqrt(2 * kappa / s^2 * (1 / qe[["q_P"]] + 2 * s * w / qe[["q_PG"]] +
                            (s * w)^2 / qe[["q_G"]]))
}

#' Standard error of a distal pair covariance estimate (eq. 9.1)
#' @param wi,wj Enrichment of the two clumps.
#' @param C Covariance \eqn{C_{ij} = \rho_{ij}\sqrt{w_i w_j}}.
#' @inheritParams se_lrcq
#' @export
se_lrcp <- function(wi, wj = wi, C = 0, n, h2, M = 1.2e6, qe) {
  qe <- qe_vec(qe)
  s <- n * h2 / M
  sqrt((1 / qe[["q_P"]] + s * (wi + wj) / qe[["q_PG"]] +
          s^2 * (wi * wj + C^2) / qe[["q_G"]]) / s^2)
}

qe_vec <- function(qe) {
  if (length(qe) == 1) qe <- c(q_P = qe, q_G = qe, q_PG = qe)
  qe
}

#' Power of per-clump LRCQ and pair-level LRCP
#'
#' \code{power_lrcq()}: one-sided test of \eqn{w = 0} at level
#' \code{alpha}. \code{power_lrcp()}: two-sided test of \eqn{\rho = 0} for two
#' clumps of equal enrichment, Bonferroni over the \eqn{K(K-1)/2} pairs among
#' \code{K} prescreened clumps.
#'
#' @param w Enrichment under the alternative.
#' @param rho Genetic-effect correlation under the alternative.
#' @param K Number of prescreened clumps.
#' @param alpha Significance level (per test for LRCQ; family-wise for LRCP).
#' @inheritParams se_lrcq
#' @return Power.
#' @export
power_lrcq <- function(w, n, h2, qe, M = 1.2e6, kappa = 1.4, alpha = 5e-7) {
  se0 <- se_lrcq(0, n, h2, M, qe, kappa)
  se1 <- se_lrcq(w, n, h2, M, qe, kappa)
  stats::pnorm((w - stats::qnorm(1 - alpha) * se0) / se1)
}

#' @rdname power_lrcq
#' @export
power_lrcp <- function(rho, w, n, h2, qe, M = 1.2e6, K = 100, alpha = 0.05) {
  a <- alpha / (K * (K - 1) / 2)
  C <- rho * w
  se0 <- se_lrcp(w, w, 0, n, h2, M, qe)
  se1 <- se_lrcp(w, w, C, n, h2, M, qe)
  zc <- stats::qnorm(1 - a / 2)
  stats::pnorm((abs(C) - zc * se0) / se1) + stats::pnorm((-abs(C) - zc * se0) / se1)
}

#' Detectable effect at a target power
#'
#' Smallest enrichment (LRCQ) or correlation (LRCP) detected with probability
#' \code{power}.
#' @param type \code{"lrcq"} (solve for w) or \code{"lrcp"} (solve for rho
#'   given \code{w}).
#' @param power Target power.
#' @param ... Passed to [power_lrcq()] or [power_lrcp()].
#' @export
detectable_effect <- function(type = c("lrcq", "lrcp"), power = 0.8, ...) {
  type <- match.arg(type)
  if (type == "lrcq") {
    f <- function(x) power_lrcq(x, ...) - power
    upper <- 1e4
  } else {
    f <- function(x) power_lrcp(x, ...) - power
    upper <- 1
  }
  if (f(upper) < 0) return(NA_real_)
  stats::uniroot(f, c(1e-8, upper))$root
}
