# Stage 3: LRCQ estimation of the genetic-effect enrichment w.
#
# Regression model (chapter 6, LRCQ key equation, R_beta = 0):
#   y_ka = z_ka^2 - C_aa = s_a (D w)_k + e_ka,   s_a = n_a h2_a / M,  D = R o R.
# Stacking traits gives the design X = s (x) D, so (fast OLS/WLS derivation)
#   X' Omega X = D' diag(omega) D,   omega_k = sum_a s_a^2 / v_ka
#   X' Omega y = D' u,               u_k     = sum_a s_a y_ka / v_ka
# and the q-trait regression collapses to one m-by-m solve.

#' Rectify negative enrichment estimates
#'
#' Chapter 6.3.2 rectification-adjustment methods, applied only when the
#' vector contains negative values.
#' \describe{
#'   \item{A}{zero every element \eqn{\le |\min(\hat w)|}.}
#'   \item{B}{zero all negatives plus the same number of smallest positives.}
#'   \item{C}{zero all elements up to the cutoff \eqn{c} at which the sum of
#'     positives in \eqn{(0, c]} is closest to the absolute sum of negatives
#'     (ties to the lower cutoff).}
#' }
#' @param w Estimated enrichment.
#' @param method One of \code{"none"}, \code{"A"}, \code{"B"}, \code{"C"},
#'   \code{"zero"} (just truncate negatives at 0).
#' @export
rectify_w <- function(w, method = c("none", "A", "B", "C", "zero")) {
  method <- match.arg(method)
  if (method == "none" || !any(w < 0, na.rm = TRUE)) return(w)
  neg <- which(w < 0)
  out <- w
  if (method == "zero") {
    out[neg] <- 0
  } else if (method == "A") {
    out[w <= abs(min(w, na.rm = TRUE))] <- 0
  } else if (method == "B") {
    pos <- which(w > 0)
    k <- min(length(neg), length(pos))
    small <- pos[order(w[pos])][seq_len(k)]
    out[c(neg, small)] <- 0
  } else if (method == "C") {
    total_neg <- -sum(w[neg])
    pos <- which(w > 0)
    sp <- sort(w[pos])
    cs <- cumsum(sp)
    # cutoff whose cumulative positive mass is closest to the negative mass
    # (ties go to the lower cutoff)
    j <- which.min(abs(cs - total_neg))
    cutoff <- sp[j]
    out[w <= cutoff] <- 0
  }
  out
}

#' LRCQ enrichment estimate within one LD window
#'
#' @param Z m-by-q Z-scores of the window SNPs.
#' @param R m-by-m LD matrix of the window.
#' @param n Length-q sample sizes.
#' @param gcov q-by-q genetic covariance (or length-q heritabilities).
#' @param M Number of SNPs heritability refers to (genome-wide).
#' @param intercept q-by-q intercept matrix (or length-q univariate
#'   intercepts). Responses are \eqn{z^2 - C_{aa}}.
#' @param method \code{"ols"}, \code{"wls"} (one reweighting step from OLS),
#'   \code{"irls"} (iterate weights to convergence), or \code{"equal"}.
#'   WLS weights are \eqn{1/\mathrm{Var}(z^2) = 1/(2\sigma^4)}. OLS/WLS
#'   estimate the power-weighted mean enrichment
#'   \eqn{\sum_a s_a^2 w_a / \sum_a s_a^2} when enrichment differs by trait;
#'   \code{"equal"} gives \eqn{\hat w = q^{-1}\sum_a D^{-1}(z_a^2 - C_{aa})/s_a},
#'   unbiased for the simple mean over traits but noisier.
#' @param ridge Ridge penalty as a fraction of the mean diagonal of
#'   \eqn{D'\Omega D}; 0 gives the exact fast OLS/WLS solution.
#' @param N_ref Reference panel size for the r-squared bias correction.
#' @param Rb Optional genetic-effect correlation for the window. When given,
#'   the chapter's cross terms are included as an offset evaluated at the
#'   current \eqn{\hat w} (\code{cross_terms = TRUE} in [lrcq()]).
#' @param pair_sum See [sigma_w()].
#' @param se \code{"model"}: sandwich SE under Gaussian Z with the full
#'   Isserlis covariance across SNPs and traits; \code{"none"}.
#' @param max_iter,tol IRLS controls.
#' @param w_plugin Enrichment used to evaluate the SE formula (per tag, or
#'   per SNP of the window); \code{NULL} uses \code{se_rectify}.
#' @param se_rectify Plug-in for the SE when \code{w_plugin} is \code{NULL}.
#'   \code{"bgA"} (default since 0.2.0, as in the LRCQ paper; theory
#'   supplement S3.10a): background-filled rule A, estimates above
#'   \eqn{|\min \hat w|} are kept and every estimate rule A would zero gets
#'   the mean of that zeroed set, floored at 0. A method of [rectify_w()]
#'   plugs in the rectified estimate instead; \code{"A"} (the 0.1 default)
#'   drops the polygenic background and makes SEs 10-25\% too small under
#'   dense enrichment. \code{"eb"}: a common window enrichment plus tag
#'   deviations shrunk by their signal-to-noise ratio. \code{"one"}: w = 1 for
#'   every SNP (simulation checks).
#' @param tags Optional indices of tag SNPs carrying the parameters
#'   (clump-level LRCQ). Responses still come from all SNPs; the estimand of
#'   each tag is its clump total \eqn{w_T + \Pi w_{\bar T}},
#'   \eqn{\Pi = (D_{\cdot T}'D_{\cdot T})^{-1}D_{\cdot T}'D_{\cdot \bar T}}
#'   (theory supplement S3.11). Use with [prune_ld()] at \eqn{r^2 < 0.5}.
#' @param psd Project \eqn{D} onto the PSD cone (clip negative eigenvalues)
#'   before solving.
#' @param skew With \code{method = "ols"}, also return the third cumulant
#'   \code{k3} of each tag's estimate for skewness-corrected p-values
#'   ([p_skew()]).
#' @param skew_tol Support truncation for the third cumulant (0 = exact).
#' @return List with \code{w} and \code{se} (one per tag), \code{iter},
#'   \code{tags}.
#' @export
lrcq_window <- function(Z, R, n, gcov, M, intercept = NULL,
                        method = c("wls", "ols", "irls", "equal"), ridge = 0,
                        N_ref = Inf, Rb = NULL, pair_sum = c("once", "chapter"),
                        se = c("model", "none"), max_iter = 20, tol = 1e-6,
                        w_plugin = NULL, tags = NULL, psd = TRUE,
                        skew = FALSE, skew_tol = 1e-3, se_rectify = "bgA") {
  method <- match.arg(method)
  if (skew && method != "ols")
    stop("skew-corrected p-values are implemented for method = \"ols\" (the LRCQ paper's estimator)")
  se <- match.arg(se)
  pair_sum <- match.arg(pair_sum)
  Z <- as.matrix(Z)
  m <- nrow(Z)
  q <- ncol(Z)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  h2 <- diag(gcov)
  s <- n * h2 / M
  Cd <- diag(intercept)
  R <- as.matrix(R)
  D <- ld_r2(R, N_ref)
  if (psd) {
    # the r2 bias correction (and float rounding) can break positive
    # semi-definiteness; project to the PSD cone before any inverse
    e <- eigen(D, symmetric = TRUE)
    if (min(e$values) < 0) D <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
  }
  if (is.null(tags)) tags <- seq_len(m)
  p <- length(tags)
  DT <- D[, tags, drop = FALSE]
  expand <- function(wt) { x <- numeric(m); x[tags] <- wt; x }
  Y <- Z^2 - matrix(Cd, m, q, byrow = TRUE)

  offset_fn <- function(w) {
    if (is.null(Rb)) return(rep(0, m))
    S <- sigma_w(pmax(expand(w), 0), Rb, pair_sum)
    diag(S) <- 0
    rowSums((R %*% S) * R)
  }
  solve_w <- function(V, off) {
    Omega <- matrix(s^2, m, q, byrow = TRUE) / V
    omega <- rowSums(Omega)
    u <- rowSums((Y - outer(off, s)) * matrix(s, m, q, byrow = TRUE) / V)
    A <- crossprod(DT, DT * omega)
    if (ridge > 0) A <- A + diag(ridge * mean(diag(A)), p)
    list(w = as.vector(solve(A, crossprod(DT, u))), A = A)
  }
  var_fn <- function(w, off) {
    mu <- outer(pmax(as.vector(DT %*% w) + off, 0), s) +
      matrix(Cd, m, q, byrow = TRUE)
    2 * mu^2
  }

  V0 <- if (method == "equal") matrix(s^2, m, q, byrow = TRUE) else matrix(1, m, q)
  off <- rep(0, m)
  fit <- solve_w(V0, off)
  iter <- 1
  if (!is.null(Rb)) {
    off <- offset_fn(fit$w)
    fit <- solve_w(V0, off)
  }
  if (method %in% c("wls", "irls")) {
    n_it <- if (method == "wls") 1 else max_iter
    for (it in seq_len(n_it)) {
      w_old <- fit$w
      V <- var_fn(fit$w, off)
      fit <- solve_w(V, off)
      if (!is.null(Rb)) off <- offset_fn(fit$w)
      iter <- iter + 1
      if (max(abs(fit$w - w_old)) < tol * max(1, max(abs(w_old)))) break
    }
  }
  w_hat <- fit$w
  # expected estimate when every SNP has w = 1 (the null of no enrichment for
  # a clump: its estimand is w_T + Pi w_Tbar, not the clump size)
  Vn <- if (method %in% c("ols", "equal")) V0 else var_fn(w_hat, off)
  om1 <- rowSums(matrix(s^2, m, q, byrow = TRUE) / Vn)
  mu1 <- as.vector(solve(fit$A, crossprod(DT, om1 * rowSums(D))))
  se_w <- rep(NA_real_, p)
  if (se == "model") {
    V <- if (method %in% c("ols", "equal")) V0 else var_fn(w_hat, off)
    Fm <- matrix(s, m, q, byrow = TRUE) / V
    cmat <- scale_matrix(n, gcov, M)
    Ainv <- solve(fit$A)
    G <- Ainv %*% t(DT)
    se_of <- function(wfull) {
      Ak <- R %*% sigma_w(pmax(wfull, 0), Rb, pair_sum) %*% R
      Vu <- 2 * (Ak^2 * (Fm %*% (cmat^2) %*% t(Fm)) +
                   2 * Ak * R * (Fm %*% (cmat * intercept) %*% t(Fm)) +
                   R^2 * (Fm %*% (intercept^2) %*% t(Fm)))
      list(Ak = Ak, se = sqrt(pmax(rowSums((G %*% Vu) * G), 0)))
    }
    if (!is.null(w_plugin)) {
      wfull <- if (length(w_plugin) == m) w_plugin else expand(w_plugin)
    } else if (se_rectify == "one") {
      wfull <- rep(1, m)  # oracle under no enrichment (simulation checks)
    } else if (se_rectify == "eb") {
      # Empirical-Bayes plug-in: a common enrichment c for every SNP of the
      # window (matched to the estimates through mu1), plus the tag
      # deviations shrunk by their signal-to-noise ratio. Clipped or
      # thresholded raw estimates bias the SE when per-SNP signal is strong.
      c0 <- max(sum(w_hat) / sum(mu1), 0)
      se0 <- se_of(rep(c0, m))$se
      d <- w_hat - c0 * mu1
      lam <- max(0, 1 - mean(se0^2) / mean(d^2))
      wfull <- rep(c0, m); wfull[tags] <- c0 + lam * d
    } else if (se_rectify == "bgA") {
      # background-filled rule A (LRCQ paper Methods; theory supplement
      # S3.10a): estimates above |min w_hat| are kept, every estimate rule A
      # would zero gets the mean of that zeroed set, floored at 0
      keep <- w_hat > abs(min(w_hat))
      bg <- if (any(!keep)) max(mean(w_hat[!keep]), 0) else 0
      wfull <- expand(ifelse(keep, w_hat, bg))
    } else wfull <- expand(rectify_w(w_hat, se_rectify))
    so <- se_of(wfull); Ak <- so$Ak; se_w <- so$se
  }
  k3 <- NULL
  if (skew && se == "model")
    k3 <- lrcq_k3(G, R, Ak, s, intercept, scale_matrix(n, gcov, M), skew_tol)
  list(w = w_hat, se = se_w, k3 = k3, mu1 = mu1, iter = iter, tags = tags)
}

# Third cumulant of the OLS LRCQ estimate per tag (theory: W_t = g_t' u with
# u_i = sum_a s_a (z_ia^2 - C_aa) a quadratic form in Gaussian Z). Both terms of
# Cov(Z) are Kronecker products, so each trace factorises into a trait-side trace
# (computed once) and a SNP-side trace over the support of g_t (entries below
# tol * max|g_t| are dropped; tol = 0 is exact).
lrcq_k3 <- function(G, R, A, s, Gam, Tm, tol = 1e-3) {
  X1 <- s * Gam; X2 <- s * Tm
  trp <- function(X, Y) sum(X * t(Y))
  X11 <- X1 %*% X1; X22 <- X2 %*% X2
  t111 <- trp(X11, X1); t112 <- trp(X11, X2); t122 <- trp(X1, X22); t222 <- trp(X22, X2)
  vapply(seq_len(nrow(G)), function(i) {
    g <- G[i, ]
    J <- if (tol > 0) which(abs(g) > tol * max(abs(g))) else seq_along(g)
    P <- g[J] * R[J, J, drop = FALSE]; Q <- g[J] * A[J, J, drop = FALSE]
    PP <- P %*% P; QQ <- Q %*% Q
    8 * (trp(PP, P) * t111 + 3 * trp(PP, Q) * t112 + 3 * trp(P, QQ) * t122 + trp(QQ, Q) * t222)
  }, 0)
}

#' Skewness-corrected upper-tail p-value
#'
#' Pearson type III (shifted gamma) approximation with mean \code{mu0},
#' variance \code{se^2} and third cumulant \code{k3}; reduces to the normal
#' reference when the skewness is negligible. The OLS LRCQ estimate is a
#' quadratic form in Gaussian Z-scores and is right-skewed when the traits have
#' few effective dimensions (sample overlap, genetic correlation), so the normal
#' reference is conservative in the bulk and anti-conservative in the upper tail.
#'
#' @param x Estimates.
#' @param mu0 Null mean (0 for "w > 0"; the clump size for "w > 1").
#' @param se Standard errors.
#' @param k3 Third cumulants (from [lrcq()] with \code{skew = TRUE}).
#' @return Upper-tail p-values.
#' @export
p_skew <- function(x, mu0, se, k3) {
  g <- k3 / se^3
  z <- (x - mu0) / se
  p <- stats::pnorm(z, lower.tail = FALSE)
  ok <- is.finite(g) & g > 1e-3
  shape <- 4 / g[ok]^2
  p[ok] <- stats::pgamma(shape + z[ok] * sqrt(shape), shape, lower.tail = FALSE)
  p
}

#' Genome-wide LRCQ: genetic-effect enrichment for every SNP
#'
#' Intercepts \eqn{C_{aa}} must come from stage 2 (genome-wide LDSC) and are
#' held fixed: with free per-trait intercepts inside a window the design
#' \eqn{[s \otimes D, I_q \otimes 1]} is one rank short and \eqn{w} is not
#' identifiable. From squared Z-scores alone only \eqn{diag(R S_w R)} is
#' identified, so \eqn{\hat w} estimates \eqn{w} under the working assumption
#' that local genetic-effect correlation is zero; otherwise its bias is
#' \eqn{D^{-1} diag(R\, C\, R)} with \eqn{C} the off-diagonal part of
#' \eqn{S_w} (theory supplement S1.5, S1.6).
#'
#' Runs [lrcq_window()] over sliding windows of LD blocks (chapter 6 sliding
#' window strategy) and keeps the core-block estimates.
#'
#' @param Z m-by-q Z-score matrix aligned to the LD reference SNPs.
#' @param ld An \code{lrcpq_ld} object.
#' @param window_size Blocks per window (odd).
#' @param rectify Rectification method for negative estimates, applied
#'   genome-wide after estimation (see [rectify_w()]).
#' @param cross_terms Include the chapter's \eqn{r_\beta} cross terms using
#'   \code{Rb} (default \code{FALSE}: chi-square-only LRCQ as in the
#'   abstract and Table 6.1).
#' @param Rb m-by-m genetic-effect correlation (needed if
#'   \code{cross_terms = TRUE}).
#' @param chr Optional chromosome per SNP (windows do not cross chromosomes).
#' @param tag_r2 If set (e.g. 0.5), estimate clump-level enrichment on a
#'   greedy-pruned tag set per window (recommended for real HapMap3 LD, where
#'   \eqn{D} is near singular); non-tag SNPs get \code{NA} and a \code{tag}
#'   column naming their tag.
#' @param windows Optional precomputed window list (\code{idx}, \code{core}),
#'   e.g. from an LD reference stored as overlapping windows.
#' @param ... Passed to [lrcq_window()] (e.g. \code{ridge}, \code{psd}).
#' @inheritParams lrcq_window
#' @return Data frame with \code{w_raw}, \code{w} (rectified), \code{se},
#'   \code{w_null1} (the expected estimate when no SNP is enriched; with
#'   tags this is \eqn{1 + \Pi 1}, about the clump size),
#'   \code{z} (\code{w_raw / se}) and \code{block}. With \code{skew = TRUE}
#'   (passed through \code{...}, OLS only) also \code{k3}, \code{skew} and the
#'   skewness-corrected one-sided p-values \code{p_gt0} (H0: w = 0) and
#'   \code{p_gt1} (H0: no enrichment, every SNP w = 1, i.e. the tag's
#'   estimand equals \code{w_null1}). Use these, not the
#'   normal reference, for real phenome-wide traits.
#' @export
lrcq <- function(Z, ld, n, gcov, M, intercept = NULL, window_size = 3,
                 method = c("wls", "ols", "irls", "equal"), rectify = "none",
                 cross_terms = FALSE, Rb = NULL, chr = NULL, tag_r2 = NULL,
                 windows = NULL, ...) {
  method <- match.arg(method)
  Z <- as.matrix(Z)
  stopifnot(nrow(Z) == ld$m)
  m <- nrow(Z)
  w <- se <- k3 <- mu1 <- rep(NA_real_, m)
  skew <- isTRUE(list(...)$skew)
  tag_of <- seq_len(m)
  if (is.null(windows)) windows <- make_windows(ld$block_id, window_size, chr)
  for (win in windows) {
    Rw <- ld$get(win$idx)
    Rbw <- if (cross_terms) as.matrix(Rb[win$idx, win$idx]) else NULL
    tags <- if (is.null(tag_r2)) NULL else prune_ld(Rw, tag_r2)
    f <- lrcq_window(Z[win$idx, , drop = FALSE], Rw, n, gcov, M, intercept,
                     method = method, N_ref = ld$N, Rb = Rbw, tags = tags, ...)
    tg <- f$tags
    core_tags <- which(tg %in% win$core)
    w[win$idx[tg[core_tags]]] <- f$w[core_tags]
    se[win$idx[tg[core_tags]]] <- f$se[core_tags]
    mu1[win$idx[tg[core_tags]]] <- f$mu1[core_tags]
    if (skew) k3[win$idx[tg[core_tags]]] <- f$k3[core_tags]
    if (!is.null(tag_r2)) {
      # each core SNP points to the tag it is most correlated with
      nearest <- tg[apply(abs(Rw[win$core, tg, drop = FALSE]), 1, which.max)]
      tag_of[win$idx[win$core]] <- win$idx[nearest]
    }
  }
  out <- data.frame(snp = seq_len(m), block = ld$block_id, w_raw = w,
                    w = NA_real_, se = se, z = w / se)
  ok <- !is.na(w)
  out$w[ok] <- rectify_w(w[ok], rectify)
  if (!is.null(tag_r2)) out$tag <- tag_of
  out$w_null1 <- mu1
  if (skew) {
    out$k3 <- k3
    out$skew <- k3 / se^3
    out$p_gt0 <- p_skew(w, 0, se, k3)
    out$p_gt1 <- p_skew(w, mu1, se, k3)
  }
  out
}
#' Brute-force stacked regression for checking the fast LRCQ solution
#'
#' Builds the \eqn{mq \times m} design \eqn{X = s \otimes D} explicitly and
#' calls \code{lm.wfit}. Only for small problems and tests.
#' @inheritParams lrcq_window
#' @param V Optional m-by-q matrix of response variances (WLS weights 1/V).
#' @export
lrcq_bruteforce <- function(Z, R, n, gcov, M, intercept = NULL, V = NULL) {
  Z <- as.matrix(Z)
  m <- nrow(Z)
  q <- ncol(Z)
  h2 <- if (is.null(dim(gcov))) gcov else diag(gcov)
  Cd <- if (is.null(intercept)) rep(1, q) else if (is.null(dim(intercept))) intercept else diag(intercept)
  s <- n * h2 / M
  D <- as.matrix(R)^2
  X <- kronecker(matrix(s, q, 1), D)
  y <- as.vector(Z^2 - matrix(Cd, m, q, byrow = TRUE))
  wt <- if (is.null(V)) rep(1, m * q) else 1 / as.vector(V)
  stats::lm.wfit(X, y, wt)$coefficients
}
