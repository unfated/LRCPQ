#' Genome-wide LRCP screen with a simulated-null calibration
#'
#' Conditional test of [lrcp_gene()] (theory S8.2) for every distal pair of
#' candidate tags, vectorised with the per-pair conditional variance and the
#' conservative orientation (the larger of the two conditional variances).
#' Each window's genetic part is the plug-in local moment ([local_moment()]),
#' estimated from the same Z-scores as the numerator, so under H0 the
#' statistic is self-normalised and light-tailed (SD below 1) and the normal
#' reference is conservative. With \code{n_null > 0} every window is redrawn
#' from the fitted separable null \eqn{Z^* \sim N(0, \Gamma\otimes R +
#' T\otimes \hat G_{A+})}, the whole screen is recomputed (local moments
#' re-estimated from \eqn{Z^*}), and the pooled null of all distal pairs gives
#' calibrated p-values (\code{p_cal}) and BH q-values (\code{q_cal}).
#'
#' @param Z Candidate tags by q Z-scores.
#' @param window Window label per tag; tags in one window share an LD matrix.
#' @param R_list Named list of tag LD matrices, one per window (names match
#'   \code{window}; rows in the order the tags appear in \code{Z}).
#' @param chr,pos Chromosome and position per tag (pairs closer than
#'   \code{min_distance} on one chromosome are not tested).
#' @param n,gcov,M,intercept As in [lrcq()]: total N, genome-wide LDSC genetic
#'   covariance (never diagonal), M = number of SNPs whose LD enters R, and
#'   fixed stage-2 intercepts.
#' @param min_distance Minimum distance (bp) for a cis pair to count as distal.
#' @param n_null Simulated null replicates of the whole screen (0 to skip).
#' @param p_report Pairs with normal p below this are returned (all pairs are
#'   still used for BH).
#' @param seed Seed for the null replicates.
#' @return List with \code{pairs} (i, j, z, p_norm, q_norm and, with
#'   \code{n_null > 0}, p_cal, q_cal for pairs with \code{p_norm <
#'   p_report}) and \code{summary} (pair count, SD of z, median lambda, and
#'   the null SD and tail fractions).
#' @export
lrcp_screen <- function(Z, window, R_list, chr, pos, n, gcov, M, intercept,
                        min_distance = 5e6, n_null = 0, p_report = 1e-3, seed = 1) {
  Z <- as.matrix(Z); Z[is.na(Z)] <- 0
  check_gcov(gcov, ncol(Z), intercept)
  q <- ncol(Z); s <- n * diag(gcov) / M; ss <- sum(s^2)
  Tm <- scale_matrix(n, gcov, M)
  wins <- unique(window)
  WD <- lapply(wins, function(w) {
    S <- which(window == w); R <- as.matrix(R_list[[as.character(w)]])
    stopifnot(nrow(R) == length(S))
    G <- local_moment(Z[S, , drop = FALSE], R, n, gcov, M, intercept)
    list(S = S, R = R, Ri = solve(R), Rh = psd_factor(R), Gh = psd_factor(G))
  })
  distal <- outer(chr, chr, "!=") | abs(outer(pos, pos, "-")) >= min_distance
  ut <- which(upper.tri(distal) & distal); rm(distal)
  if (!length(ut)) stop("no distal pairs")
  screen_z <- function(Zx) {
    m <- nrow(Zx); X <- matrix(0, m, q); tauR <- MG <- numeric(m)
    for (w in WD) {
      Zs <- Zx[w$S, , drop = FALSE]
      G <- local_moment(Zs, w$R, n, gcov, M, intercept)
      X[w$S, ] <- w$Ri %*% Zs; tauR[w$S] <- diag(w$Ri)
      MG[w$S] <- diag(w$Ri %*% G %*% w$Ri)
    }
    sY <- sweep(X, 2, s, "*")
    KG <- rowSums((sY %*% intercept) * sY); KT <- rowSums((sY %*% Tm) * sY)
    VA <- (outer(tauR, KG) + outer(MG, KT)) / ss^2
    (tcrossprod(X, sY) / ss)[ut] / sqrt(pmax(VA, t(VA))[ut])
  }
  z <- screen_z(Z)
  p_norm <- 2 * stats::pnorm(-abs(z))
  ij <- arrayInd(ut, rep(nrow(Z), 2))
  summ <- data.frame(n_pairs = length(z), sd_z = stats::sd(z),
                     lambda_median = stats::median(z^2) / stats::qchisq(0.5, 1),
                     frac_p05 = mean(p_norm < 0.05), frac_p1e3 = mean(p_norm < 1e-3))
  keep <- which(p_norm < p_report)
  pairs <- data.frame(i = ij[keep, 1], j = ij[keep, 2], z = z[keep], p_norm = p_norm[keep],
                      q_norm = stats::p.adjust(p_norm, "BH")[keep])
  if (n_null > 0) {
    Ch <- psd_factor(intercept); Th <- psd_factor(Tm)
    brk <- seq(0, 12, by = 0.001); h0 <- numeric(length(brk)); tail0 <- numeric(0)
    nsd <- numeric(n_null)
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    for (b in seq_len(n_null)) {
      Zs <- matrix(0, nrow(Z), q)
      for (w in WD) {
        m1 <- length(w$S)
        Zs[w$S, ] <- w$Rh %*% matrix(stats::rnorm(m1 * q), m1) %*% t(Ch) +
          w$Gh %*% matrix(stats::rnorm(m1 * q), m1) %*% t(Th)
      }
      z0 <- abs(screen_z(Zs)); nsd[b] <- sqrt(mean(z0^2))
      h0 <- h0 + tabulate(pmin(findInterval(z0, brk), length(brk)), nbins = length(brk))
      tail0 <- c(tail0, z0[z0 > 4])
    }
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
    p_cal <- calibrate_p(abs(z), brk, h0, tail0)
    pairs$p_cal <- p_cal[keep]; pairs$q_cal <- stats::p.adjust(p_cal, "BH")[keep]
    summ$null_sd <- mean(nsd); summ$null_reps <- n_null
    summ$n_q05_norm <- sum(stats::p.adjust(p_norm, "BH") < 0.05)
    summ$n_q05_cal <- sum(stats::p.adjust(p_cal, "BH") < 0.05)
  }
  list(pairs = pairs[order(pairs$p_norm), ], summary = summ)
}

# Empirical P(|z*| >= a) from a binned null (counts h on breaks brk) with the
# exact null values above 4 kept in tail; (1 + count) / (N + 1).
calibrate_p <- function(a, brk, h, tail) {
  N <- sum(h); cum <- rev(cumsum(rev(h)))
  b <- pmin(findInterval(a, brk), length(brk))
  out <- (1 + cum[b]) / (N + 1)
  k <- a > 4
  if (any(k)) {
    r <- sort(tail)
    out[k] <- (1 + length(r) - findInterval(a[k], r, left.open = TRUE)) / (N + 1)
  }
  out
}
