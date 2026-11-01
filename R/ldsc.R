# Stage 1 and 2: LD score regression for heritability, co-heritability and
# intercepts (Bulik-Sullivan et al. 2015a, 2015b).

#' Weighted simple regression with block jackknife
#'
#' Fits \eqn{y = a + b x} by weighted least squares and returns delete-one-block
#' jackknife standard errors computed from block sufficient statistics.
#'
#' @param x,y Predictor and response.
#' @param wt Regression weights.
#' @param block Jackknife block label per observation.
#' @param intercept Fixed intercept value, or \code{NULL} to estimate it.
#' @return List with \code{coef} (intercept, slope), \code{se}, and the
#'   jackknife pseudo-values matrix \code{jk}.
#' @keywords internal
wls_jackknife <- function(x, y, wt, block, intercept = NULL) {
  fixed <- !is.null(intercept)
  if (fixed) y <- y - intercept
  X <- if (fixed) cbind(x) else cbind(1, x)
  p <- ncol(X)
  f <- factor(block)
  nb <- nlevels(f)
  # per-block X'WX and X'Wy
  XtWX_b <- array(0, c(p, p, nb))
  XtWy_b <- matrix(0, p, nb)
  for (j in seq_len(p)) {
    XtWy_b[j, ] <- rowsum(wt * X[, j] * y, f, reorder = TRUE)[, 1]
    for (k in j:p) {
      v <- rowsum(wt * X[, j] * X[, k], f, reorder = TRUE)[, 1]
      XtWX_b[j, k, ] <- v
      XtWX_b[k, j, ] <- v
    }
  }
  XtWX <- apply(XtWX_b, c(1, 2), sum)
  XtWy <- rowSums(XtWy_b)
  est <- solve(XtWX, XtWy)
  jk <- t(vapply(seq_len(nb), function(b) {
    solve(XtWX - XtWX_b[, , b], XtWy - XtWy_b[, b])
  }, numeric(p)))
  if (p == 1) jk <- matrix(jk, ncol = 1)
  pseudo <- nb * matrix(est, nb, p, byrow = TRUE) - (nb - 1) * jk
  se <- sqrt(apply(pseudo, 2, stats::var) / nb)
  if (fixed) {
    coef <- c(intercept, est)
    se <- c(0, se)
  } else {
    coef <- est
  }
  names(coef) <- names(se) <- c("intercept", "slope")
  list(coef = coef, se = se, jk = jk)
}

#' Univariate LD score regression
#'
#' \eqn{E[\chi^2_k] = N h^2 l_k / M + a}. Weights follow LDSC: the inverse of
#' the regression-SNP LD score times the heteroscedastic variance
#' \eqn{2(N h^2 l/M + a)^2}, updated once from a first fit.
#'
#' @param z Z-scores (length m).
#' @param ldscore LD scores (length m).
#' @param n Sample size (scalar or per SNP).
#' @param M Number of SNPs heritability refers to.
#' @param w_ld Regression-SNP LD scores for weighting (default \code{ldscore}).
#' @param intercept Fixed intercept, or \code{NULL} to estimate.
#' @param n_blocks Number of contiguous jackknife blocks.
#' @param chisq_max Drop SNPs with larger chi-square (LDSC default:
#'   \code{max(80, 0.001 * N)}).
#' @return List with \code{h2}, \code{h2_se}, \code{intercept},
#'   \code{intercept_se}.
#' @export
ldsc_h2 <- function(z, ldscore, n, M, w_ld = ldscore, intercept = NULL,
                    n_blocks = 200, chisq_max = NULL) {
  m <- length(z)
  n <- rep_len(n, m)
  chi2 <- z^2
  if (is.null(chisq_max)) chisq_max <- max(80, 0.001 * max(n))
  ok <- is.finite(chi2) & chi2 <= chisq_max
  x <- n * ldscore / M
  block <- ceiling(seq_len(m) / ceiling(m / n_blocks))
  wt_fun <- function(h2, a) {
    h2 <- min(max(h2, 0), 1)
    1 / (pmax(w_ld, 1) * 2 * (a + h2 * x)^2)
  }
  h2_0 <- (mean(chi2[ok]) - 1) / mean(x[ok])
  a0 <- if (is.null(intercept)) 1 else intercept
  fit <- wls_jackknife(x[ok], chi2[ok], wt_fun(h2_0, a0)[ok], block[ok], intercept)
  fit <- wls_jackknife(x[ok], chi2[ok],
                       wt_fun(fit$coef[2], fit$coef[1])[ok], block[ok], intercept)
  list(h2 = unname(fit$coef[2]), h2_se = unname(fit$se[2]),
       intercept = unname(fit$coef[1]), intercept_se = unname(fit$se[1]))
}

#' Cross-trait LD score regression
#'
#' \eqn{E[z_{ka} z_{kb}] = \sqrt{n_a n_b} h_{ab} l_k/M + \rho_{ab}}, where the
#' intercept \eqn{\rho_{ab} = r_{ab} o_{ab}} is the stage-2 parameter.
#'
#' @param z1,z2 Z-scores of the two traits.
#' @param n1,n2 Sample sizes.
#' @param h2_1,h2_2,int_1,int_2 Univariate fits used for the weights.
#' @inheritParams ldsc_h2
#' @return List with \code{gcov} (co-heritability), \code{gcov_se},
#'   \code{intercept}, \code{intercept_se}, and the genetic correlation
#'   \code{rg}.
#' @export
ldsc_gcov <- function(z1, z2, ldscore, n1, n2, M, h2_1, h2_2, int_1 = 1,
                      int_2 = 1, w_ld = ldscore, intercept = NULL,
                      n_blocks = 200) {
  m <- length(z1)
  n1 <- rep_len(n1, m)
  n2 <- rep_len(n2, m)
  y <- z1 * z2
  ok <- is.finite(y)
  x <- sqrt(n1 * n2) * ldscore / M
  block <- ceiling(seq_len(m) / ceiling(m / n_blocks))
  v1 <- n1 * max(h2_1, 0) * ldscore / M + int_1
  v2 <- n2 * max(h2_2, 0) * ldscore / M + int_2
  wt_fun <- function(g, a) 1 / (pmax(w_ld, 1) * (v1 * v2 + (a + g * x)^2))
  a0 <- if (is.null(intercept)) 0 else intercept
  fit <- wls_jackknife(x[ok], y[ok], wt_fun(0, a0)[ok], block[ok], intercept)
  fit <- wls_jackknife(x[ok], y[ok], wt_fun(fit$coef[2], fit$coef[1])[ok],
                       block[ok], intercept)
  g <- unname(fit$coef[2])
  list(gcov = g, gcov_se = unname(fit$se[2]),
       intercept = unname(fit$coef[1]), intercept_se = unname(fit$se[1]),
       rg = g / sqrt(max(h2_1, 1e-12) * max(h2_2, 1e-12)))
}

#' Stages 1 and 2 for all traits: genetic covariance and intercept matrices
#'
#' Runs [ldsc_h2()] on every trait and [ldsc_gcov()] on every pair, returning
#' the q-by-q genetic covariance matrix (stage 1) and intercept matrix
#' (stage 2) consumed by the LRCQ/LRCP estimators.
#'
#' @param Z m-by-q Z-score matrix.
#' @param n Length-q sample sizes.
#' @param pairs Estimate the off-diagonal (co-heritability, cross-trait
#'   intercept); if \code{FALSE} they are set to 0.
#' @param fix_intercept Fix univariate intercepts at 1 (simulation without
#'   stratification) instead of estimating them.
#' @inheritParams ldsc_h2
#' @return List with \code{gcov}, \code{gcov_se}, \code{intercept},
#'   \code{intercept_se} (q-by-q) and \code{h2}.
#' @export
ldsc_matrix <- function(Z, ldscore, n, M, w_ld = ldscore, pairs = TRUE,
                        fix_intercept = FALSE, n_blocks = 200) {
  q <- ncol(Z)
  G <- Gse <- C <- Cse <- matrix(0, q, q)
  for (a in seq_len(q)) {
    f <- ldsc_h2(Z[, a], ldscore, n[a], M, w_ld,
                 intercept = if (fix_intercept) 1 else NULL, n_blocks = n_blocks)
    G[a, a] <- f$h2
    Gse[a, a] <- f$h2_se
    C[a, a] <- f$intercept
    Cse[a, a] <- f$intercept_se
  }
  if (pairs && q > 1) {
    for (a in 1:(q - 1)) for (b in (a + 1):q) {
      f <- ldsc_gcov(Z[, a], Z[, b], ldscore, n[a], n[b], M, G[a, a], G[b, b],
                     C[a, a], C[b, b], w_ld, n_blocks = n_blocks)
      G[a, b] <- G[b, a] <- f$gcov
      Gse[a, b] <- Gse[b, a] <- f$gcov_se
      C[a, b] <- C[b, a] <- f$intercept
      Cse[a, b] <- Cse[b, a] <- f$intercept_se
    }
  }
  list(gcov = G, gcov_se = Gse, intercept = C, intercept_se = Cse, h2 = diag(G))
}
