# Stage 4: LRCP estimation of genetic-effect correlation R_beta.
#
# Distal regions (no LD between region 1 and region 2):
#   E[z1_a z2_b'] = f c_ab A1 P A2',  A_r = R_r W_r^{1/2}[, screened],
# with P the cross-region block of R_beta (screened SNPs), f = 1 ("once") or
# 2 ("chapter"). Stacking trait pairs with weights K gives the averaged response
#   Ybar = Z1 K Z2' / ||K||_F^2   and   P_hat = (A1'A1)^-1 A1' Ybar A2 (A2'A2)^-1 / f,
# the cross-product analogue of the fast LRCQ OLS solution.

#' Trait-pair weighting matrix for LRCP
#' @keywords internal
lrcp_K <- function(cmat, trait_pairs) {
  if (trait_pairs == "same") diag(diag(cmat), nrow(cmat)) else cmat
}

#' Screen SNPs for LRCP on estimated enrichment
#'
#' Chapter 6 prescreening: \eqn{\rho_{ij}} is set to 0 and not estimated if
#' \eqn{w_i < t} or \eqn{w_j < t}.
#' @param w Enrichment estimates.
#' @param threshold Cutoff \eqn{t}.
#' @param z Optional enrichment z-scores; if given, SNPs must also have
#'   \code{z >= z_min}.
#' @param z_min Minimum z-score.
#' @export
screen_snps <- function(w, threshold = 0.5, z = NULL, z_min = -Inf) {
  keep <- is.finite(w) & w >= threshold
  if (!is.null(z)) keep <- keep & is.finite(z) & z >= z_min
  which(keep)
}

#' LRCP between two distal regions (regression)
#'
#' Estimates the cross-region block of \eqn{R_\beta} between screened SNPs
#' of two regions that share no LD (e.g. different chromosomes, such as the
#' TP53 and MDM2 loci), regressing Z-score cross-products on LD cross
#' components.
#'
#' @param Z1,Z2 Z-scores (m1-by-q, m2-by-q) of the two regions.
#' @param R1,R2 LD matrices of the two regions.
#' @param w1,w2 Enrichment of all SNPs in each region (stage 3 output,
#'   rectified, or true values in a simulation).
#' @param S1,S2 Indices of screened SNPs in each region (default
#'   \code{w > 0}).
#' @param n,gcov,M,intercept As in [lrcq_window()].
#' @param trait_pairs Use only same-trait products \eqn{z_{ka} z_{la}}
#'   (\code{"same"}) or all trait pairs weighted by \eqn{c_{ab}}
#'   (\code{"all"}).
#' @param method \code{"gls"} (default): generalised least squares under the
#'   null covariance \eqn{R_2 \otimes R_1}, which collapses to
#'   \eqn{W^{-1/2} R_{PP}^{-1} \bar Y_{P_1 P_2} R_{PP}^{-1} W^{-1/2}} and needs
#'   only the screened SNPs' Z-scores and LD (theory supplement S1.5);
#'   \code{"ols"}, or \code{"wls"} (separable weights
#'   \eqn{1/(\sigma^2_{ka}\sigma^2_{lb})}, same-trait products only). Both
#'   report sandwich SEs that use the full q-by-q Isserlis covariance of the
#'   cross-products under \eqn{\rho = 0}, so genetic correlation and sample
#'   overlap between GWAS are accounted for. The WLS weights ignore LD between
#'   SNPs, and in simulation WLS is not more efficient than OLS; use
#'   [lrcp_mle()] for efficiency.
#' @param ridge Ridge penalty as a fraction of the mean diagonal of the
#'   normal-equation matrix.
#' @param pair_sum See [sigma_w()].
#' @param N_ref Reference panel size (r-squared bias correction for the
#'   variance terms).
#' @param denominators Scale for \eqn{\rho}: \code{"w"} uses the supplied
#'   enrichment; \code{"R"} re-estimates the candidates' variances with the
#'   same R-projection as the numerator ([local_variance()]), which keeps
#'   \eqn{|\hat\rho| \le 1} meaningful when candidates are tags in tight LD
#'   (theory supplement S4.7).
#' @param se_type \code{"null"} (default): SEs under \eqn{\rho = 0}, the
#'   right reference for testing. \code{"plugin"} (GLS and OLS): adds the
#'   Isserlis term from the cross-region covariance at \eqn{\hat\rho},
#'   \eqn{tr(KcKc)/\|K\|^4 \cdot \hat\rho^2} for unit-scaled maps, which the
#'   null SE omits; use it for confidence intervals when \eqn{|\rho|} is large
#'   (in simulation the null SE understates the SD by about 25\% at
#'   \eqn{\rho = 0.9} with 30 traits and is accurate for
#'   \eqn{|\rho| \le 0.4}).
#' @return List with \code{rho} (p1-by-p2 estimates), \code{se}, \code{z},
#'   and the screened indices.
#' @export
lrcp_distal <- function(Z1, Z2, R1, R2, w1, w2, n, gcov, M, intercept = NULL,
                        S1 = which(w1 > 0), S2 = which(w2 > 0),
                        trait_pairs = c("same", "all"),
                        method = c("gls", "ols", "wls"), ridge = 0,
                        pair_sum = c("once", "chapter"), N_ref = Inf,
                        denominators = c("w", "R"), se_type = c("null", "plugin")) {
  trait_pairs <- match.arg(trait_pairs)
  se_type <- match.arg(se_type)
  method <- match.arg(method)
  pair_sum <- match.arg(pair_sum)
  denominators <- match.arg(denominators)
  Z1 <- as.matrix(Z1); Z2 <- as.matrix(Z2)
  R1 <- as.matrix(R1); R2 <- as.matrix(R2)
  q <- ncol(Z1)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  if (denominators == "R") {
    w1[S1] <- local_variance(Z1, R1, S1, n, gcov, M, intercept)
    w2[S2] <- local_variance(Z2, R2, S2, n, gcov, M, intercept)
  }
  # candidates with a non-positive variance have no defined rho: drop them
  # (re-projecting the rest), return NA in their rows/columns, and warn
  ok1 <- w1[S1] > 0; ok2 <- w2[S2] > 0
  if (!all(ok1) || !all(ok2)) {
    warning(sum(!ok1) + sum(!ok2), " candidate(s) with non-positive ",
            if (denominators == "R") "R-projected " else "", "variance dropped; ",
            "their rho is NA", call. = FALSE)
    p1 <- length(S1); p2 <- length(S2)
    out <- list(rho = matrix(NA_real_, p1, p2), se = matrix(NA_real_, p1, p2))
    if (any(ok1) && any(ok2)) {
      f <- lrcp_distal(Z1, Z2, R1, R2, w1, w2, n, gcov, M, intercept,
                       S1 = S1[ok1], S2 = S2[ok2], trait_pairs = trait_pairs,
                       method = method, ridge = ridge, pair_sum = pair_sum,
                       N_ref = N_ref, denominators = denominators, se_type = se_type)
      i1 <- match(f$S1, S1); i2 <- match(f$S2, S2)
      out$rho[i1, i2] <- f$rho; out$se[i1, i2] <- f$se
    }
    return(list(rho = out$rho, se = out$se, z = out$rho / out$se, S1 = S1, S2 = S2))
  }
  f <- if (pair_sum == "once") 1 else 2
  w1 <- pmax(w1, 0); w2 <- pmax(w2, 0)
  A1 <- R1[, S1, drop = FALSE] %*% diag(sqrt(w1[S1]), length(S1))
  A2 <- R2[, S2, drop = FALSE] %*% diag(sqrt(w2[S2]), length(S2))
  cmat <- scale_matrix(n, gcov, M)
  K <- lrcp_K(cmat, trait_pairs)
  p1 <- length(S1); p2 <- length(S2)
  G1 <- R1 %*% (w1 * R1)   # R W R, local rho = 0
  G2 <- R2 %*% (w2 * R2)
  addridge <- function(A) if (ridge > 0) A + diag(ridge * mean(diag(A)), nrow(A)) else A

  if (method %in% c("ols", "gls")) {
    Ybar <- Z1 %*% K %*% t(Z2) / sum(K^2)
    if (method == "ols") {
      L1 <- solve(addridge(crossprod(A1)), t(A1)) / f
      L2 <- solve(addridge(crossprod(A2)), t(A2))
    } else {
      # GLS under the null Kronecker covariance R2 (x) R1 collapses to
      # W^{-1/2} R[P,P]^{-1} Ybar[P1,P2] R[P,P]^{-1} W^{-1/2}
      sel <- function(m, S) diag(m)[S, , drop = FALSE]
      L1 <- (1 / sqrt(w1[S1])) * solve(addridge(R1[S1, S1, drop = FALSE]), sel(nrow(R1), S1)) / f
      L2 <- (1 / sqrt(w2[S2])) * solve(addridge(R2[S2, S2, drop = FALSE]), sel(nrow(R2), S2))
    }
    rho <- L1 %*% Ybar %*% t(L2)
    # Isserlis covariance of vec(Ybar) under rho = 0:
    # sum over (Y, X) in {(c, G), (C, R)} of tr(Y K X K) (X2 (x) Y1) / ||K||^4
    trm <- function(Y, X) sum(diag(Y %*% K %*% X %*% K))
    V1 <- list(c = L1 %*% G1 %*% t(L1), C = L1 %*% R1 %*% t(L1))
    V2 <- list(c = L2 %*% G2 %*% t(L2), C = L2 %*% R2 %*% t(L2))
    mats <- list(c = cmat, C = intercept)
    var <- matrix(0, p1, p2)
    for (y in c("c", "C")) for (x in c("c", "C")) {
      var <- var + trm(mats[[y]], mats[[x]]) * outer(diag(V1[[y]]), diag(V2[[x]]))
    }
    if (se_type == "plugin") {
      # second Isserlis term, Cov(z1a, z2d) Cov(z2b, z1c), from the cross-region
      # moment c_ad R1 W1^1/2 Rb12 W2^1/2 R2 at the estimate
      X12 <- A1 %*% (f * rho) %*% t(A2)
      var <- var + trm(cmat, cmat) * (L1 %*% X12 %*% t(L2))^2
    }
    se <- sqrt(pmax(var, 0)) / sum(K^2)
  } else {
    if (se_type == "plugin") warning("se_type = \"plugin\" is implemented for GLS and OLS; WLS reports the null SE")
    # separable WLS over same-trait products
    s <- diag(cmat)
    v1 <- outer(diag(G1), s) + matrix(diag(intercept), nrow(Z1), q, byrow = TRUE)
    v2 <- outer(diag(G2), s) + matrix(diag(intercept), nrow(Z2), q, byrow = TRUE)
    LHS <- matrix(0, p1 * p2, p1 * p2)
    rhs <- matrix(0, p1, p2)
    for (a in seq_len(q)) {
      K1 <- crossprod(A1, A1 / v1[, a])
      K2 <- crossprod(A2, A2 / v2[, a])
      LHS <- LHS + s[a]^2 * kronecker(K2, K1)
      rhs <- rhs + s[a] * crossprod(A1 / v1[, a], tcrossprod(Z1[, a], Z2[, a])) %*% (A2 / v2[, a])
    }
    LHS <- addridge(LHS)
    Linv <- solve(LHS)
    rho <- matrix(Linv %*% as.vector(rhs), p1, p2) / f
    # sandwich covariance under rho = 0 with the full q-by-q Isserlis blocks:
    # Cov(rhs) = sum_ab s_a s_b (B2a' S2_ab B2b) (x) (B1a' S1_ab B1b),
    # S_r,ab = c_ab G_r + C_ab R_r
    B1 <- lapply(seq_len(q), function(a) A1 / v1[, a])
    B2 <- lapply(seq_len(q), function(a) A2 / v2[, a])
    meat <- matrix(0, p1 * p2, p1 * p2)
    for (a in seq_len(q)) for (b in seq_len(q)) {
      if (cmat[a, b] == 0 && intercept[a, b] == 0) next
      M1 <- crossprod(B1[[a]], (cmat[a, b] * G1 + intercept[a, b] * R1) %*% B1[[b]])
      M2 <- crossprod(B2[[a]], (cmat[a, b] * G2 + intercept[a, b] * R2) %*% B2[[b]])
      meat <- meat + s[a] * s[b] * kronecker(M2, M1)
    }
    se <- matrix(sqrt(pmax(rowSums((Linv %*% meat) * Linv), 0)), p1, p2) / f
  }
  list(rho = rho, se = se, z = rho / se, S1 = S1, S2 = S2)
}

#' LRCP by Gaussian maximum likelihood (matrix-normal model)
#'
#' Under the matrix-normal model the stacked Z-scores of two regions satisfy
#' \eqn{Cov(vec Z) = c \otimes (R S_w R) + C \otimes R}. Simultaneously
#' diagonalising the q-by-q matrices \eqn{c} and \eqn{C} turns this into q
#' independent Gaussian vectors with covariance
#' \eqn{\lambda_t R S_w R + R}, so the exact log-likelihood costs q Cholesky
#' factorisations. The cross-region correlations \eqn{\rho_{ij}} of screened
#' SNPs are estimated by BFGS with an analytic gradient over a parameterisation
#' that keeps \eqn{R_\beta} positive semi-definite (spectral norm of the
#' cross block below 1); \code{w} is held at its stage-3 value and
#' within-region correlation at 0.
#'
#' @inheritParams lrcp_distal
#' @param start Starting values (p1-by-p2), default the OLS estimate.
#' @param bound Bound on the spectral norm of the cross-region block used to
#'   rescale an infeasible starting value.
#' @return List with \code{rho}, \code{se} (inverse observed information from
#'   a numerical Hessian of the analytic gradient), \code{z}, \code{loglik},
#'   \code{loglik0} (at \eqn{\rho = 0}) and the likelihood-ratio statistic
#'   \code{lrt} for \eqn{\rho = 0}.
#' @export
lrcp_mle <- function(Z1, Z2, R1, R2, w1, w2, n, gcov, M, intercept = NULL,
                     S1 = which(w1 > 0), S2 = which(w2 > 0),
                     pair_sum = c("once", "chapter"), start = NULL,
                     bound = 0.99) {
  pair_sum <- match.arg(pair_sum)
  Z1 <- as.matrix(Z1); Z2 <- as.matrix(Z2)
  q <- ncol(Z1)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  f <- if (pair_sum == "once") 1 else 2
  m1 <- nrow(Z1); m2 <- nrow(Z2)
  R <- as.matrix(Matrix::bdiag(as.matrix(R1), as.matrix(R2)))
  w <- pmax(c(w1, w2), 0)
  I1 <- S1; I2 <- m1 + S2
  p1 <- length(S1); p2 <- length(S2)
  cmat <- scale_matrix(n, gcov, M)
  # simultaneous diagonalisation: T C T' = I, T c T' = diag(lambda)
  eC <- eigen(intercept, symmetric = TRUE)
  Ch <- eC$vectors %*% diag(1 / sqrt(eC$values), q) %*% t(eC$vectors)
  e2 <- eigen(Ch %*% cmat %*% Ch, symmetric = TRUE)
  Tm <- t(e2$vectors) %*% Ch
  lambda <- e2$values
  Zt <- rbind(Z1, Z2) %*% t(Tm)
  sw <- sqrt(w)
  build_S <- function(P) {
    S <- diag(w, length(w))
    S[I1, I2] <- f * P * outer(sw[I1], sw[I2])
    S[I2, I1] <- t(S[I1, I2])
    S
  }
  RW <- function(P) R %*% build_S(P) %*% R
  nll_grad <- function(par, grad = TRUE) {
    P <- matrix(par, p1, p2)
    G <- RW(P)
    nll <- 0
    gr <- matrix(0, length(w), length(w))
    for (t in seq_len(q)) {
      Sig <- lambda[t] * G + R
      L <- tryCatch(chol(Sig), error = function(e) NULL)
      if (is.null(L)) return(if (grad) list(nll = 1e10, gr = rep(0, length(par))) else 1e10)
      a <- backsolve(L, forwardsolve(t(L), Zt[, t]))
      nll <- nll + sum(log(diag(L))) + 0.5 * sum(Zt[, t] * a)
      if (grad) {
        Si <- chol2inv(L)
        gr <- gr + lambda[t] * (Si - tcrossprod(a))
      }
    }
    if (!grad) return(nll)
    # d nll / d P_ij = f sqrt(w_i w_j) [R (Si - a a') R]_{ij} summed (x2 for symmetry, x1/2)
    RgR <- R %*% gr %*% R
    g <- f * RgR[I1, I2] * outer(sw[I1], sw[I2])
    list(nll = nll, gr = as.vector(g))
  }
  if (is.null(start)) {
    start <- lrcp_distal(Z1, Z2, R1, R2, w1, w2, n, gcov, M, intercept, S1, S2,
                         trait_pairs = "all", pair_sum = pair_sum)$rho
  }
  # Positive-definiteness: with within-region R_beta = I, the screened block
  # [[I, P], [P', I]] is PSD iff ||P||_2 <= 1. Optimise over unconstrained K
  # with P = K (I + K'K)^{-1/2}, a bijection onto the open unit ball.
  msqrt_inv <- function(S) {
    e <- eigen(S, symmetric = TRUE)
    e$vectors %*% (t(e$vectors) / sqrt(pmax(e$values, 1e-12)))
  }
  to_P <- function(k) {
    K <- matrix(k, p1, p2)
    as.vector(K %*% msqrt_inv(diag(p2) + crossprod(K)))
  }
  to_K <- function(P) {
    nrm <- if (length(P)) max(svd(P)$d) else 0
    if (nrm >= bound) P <- P * bound / nrm * 0.9
    as.vector(P %*% msqrt_inv(diag(p2) - crossprod(P)))
  }
  fn <- function(k) nll_grad(to_P(k), FALSE)
  gr <- function(k) {
    P <- to_P(k)
    g <- nll_grad(P, TRUE)$gr
    # chain rule through the small matrix map by central differences
    J <- vapply(seq_along(k), function(j) {
      e <- numeric(length(k)); e[j] <- 1e-6
      (to_P(k + e) - to_P(k - e)) / 2e-6
    }, numeric(length(k)))
    as.vector(crossprod(J, g))
  }
  k0 <- to_K(matrix(start, p1, p2))
  if (fn(k0) >= 1e10) k0 <- numeric(p1 * p2)
  opt <- stats::optim(k0, fn, gr, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-12))
  opt$par <- to_P(opt$par)
  par <- opt$par
  # numerical Hessian from the analytic gradient
  h <- 1e-4
  H <- vapply(seq_along(par), function(j) {
    e <- numeric(length(par)); e[j] <- h
    (nll_grad(par + e)$gr - nll_grad(par - e)$gr) / (2 * h)
  }, numeric(length(par)))
  H <- (H + t(H)) / 2
  V <- tryCatch(solve(H), error = function(e) MASS_ginv(H))
  se <- matrix(sqrt(pmax(diag(V), 0)), p1, p2)
  rho <- matrix(par, p1, p2)
  nll0 <- nll_grad(rep(0, length(par)), FALSE)
  list(rho = rho, se = se, z = rho / se, loglik = -opt$value, loglik0 = -nll0,
       lrt = 2 * (nll0 - opt$value), convergence = opt$convergence,
       S1 = S1, S2 = S2)
}

#' Moore-Penrose inverse (fallback for singular Hessians)
#' @keywords internal
MASS_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X)
  pos <- s$d > max(tol * s$d[1L], 0)
  s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

#' LRCP within one LD region (local pairs)
#'
#' Estimates \eqn{\rho_{ij}} for pairs of screened SNPs inside one LD window,
#' where LD confounding is strongest. The known diagonal part
#' \eqn{c_{ab} R W R + C_{ab} R} is removed, and the averaged residual
#' cross-product matrix is regressed on the symmetric LD cross components
#' \eqn{f (a_i a_j' + a_j a_i')}, \eqn{a_i = \sqrt{w_i} R_{\cdot i}}.
#' Intended for small screened sets (the design has p(p-1)/2 columns).
#'
#' @param Z m-by-q Z-scores of the region.
#' @param R LD matrix.
#' @param w Enrichment of the region's SNPs.
#' @param S Screened SNP indices.
#' @inheritParams lrcp_distal
#' @return Data frame of pairs (\code{i}, \code{j}, \code{rho}).
#' @export
lrcp_local <- function(Z, R, w, n, gcov, M, intercept = NULL, S = which(w > 0),
                       trait_pairs = c("same", "all"),
                       pair_sum = c("once", "chapter"), ridge = 0) {
  trait_pairs <- match.arg(trait_pairs)
  pair_sum <- match.arg(pair_sum)
  Z <- as.matrix(Z); R <- as.matrix(R)
  q <- ncol(Z); m <- nrow(Z)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  f <- if (pair_sum == "once") 1 else 2
  w <- pmax(w, 0)
  cmat <- scale_matrix(n, gcov, M)
  K <- lrcp_K(cmat, trait_pairs)
  kk <- sum(K^2)
  Ybar <- Z %*% K %*% t(Z) / kk -
    (sum(K * cmat) / kk) * (R %*% (w * R)) - (sum(K * intercept) / kk) * R
  pr <- utils::combn(S, 2)
  A <- R[, S, drop = FALSE] %*% diag(sqrt(w[S]), length(S))
  col <- function(u, v) {
    ai <- A[, match(u, S)]; aj <- A[, match(v, S)]
    as.vector(tcrossprod(ai, aj) + tcrossprod(aj, ai))
  }
  X <- f * vapply(seq_len(ncol(pr)), function(k) col(pr[1, k], pr[2, k]), numeric(m * m))
  XtX <- crossprod(X)
  if (ridge > 0) XtX <- XtX + diag(ridge * mean(diag(XtX)), ncol(X))
  rho <- solve(XtX, crossprod(X, as.vector(Ybar)))
  data.frame(i = pr[1, ], j = pr[2, ], rho = as.vector(rho))
}

#' R-projected local variance of candidate SNPs
#'
#' \eqn{\hat V_T = R_{TT}^{-1}[\sum_a s_a (z_{Ta} z_{Ta}' - C_{aa} R_{TT})]
#' R_{TT}^{-1} / \|s\|^2}; its diagonal estimates the candidates'
#' enrichment on the same projection used by the GLS numerator of
#' [lrcp_distal()].
#' @param Z,R Region Z-scores and LD.
#' @param S Candidate indices.
#' @inheritParams lrcp_distal
#' @return Length-|S| vector.
#' @export
local_variance <- function(Z, R, S, n, gcov, M, intercept = NULL) {
  Z <- as.matrix(Z)
  q <- ncol(Z)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  s <- n * diag(gcov) / M
  RT <- as.matrix(R)[S, S, drop = FALSE]
  ZT <- Z[S, , drop = FALSE]
  Ybar <- (ZT %*% (s * t(ZT)) - sum(s * diag(intercept)) * RT) / sum(s^2)
  Ri <- solve(RT)
  diag(Ri %*% Ybar %*% Ri)
}
