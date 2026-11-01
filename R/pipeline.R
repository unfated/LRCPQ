# End-to-end pipelines for the four-stage LRCQ/LRCP estimator.

#' Run LRCQ end to end (stages 1-3)
#'
#' Stage 1-2 values are taken from \code{gcov}/\code{intercept} when given
#' (e.g. genome-wide LDSC, the Pan-UKB h2 manifest, or the truth in a
#' simulation); otherwise they are estimated with [ldsc_matrix()] from
#' \code{Z} and \code{ldscore}, which needs genome-scale data. Stage 3 then
#' estimates enrichment per tag clump (default) or per SNP.
#'
#' @param Z m-by-q Z-scores aligned to \code{ld}.
#' @param ld An \code{lrcpq_ld} object.
#' @param n Sample sizes.
#' @param M Number of SNPs heritability refers to (the size of the SNP
#'   universe behind \code{D}, e.g. 1,094,844 HapMap3 SNPs; not Pan-UKB's
#'   M_5_50, see [check_scale()]). \code{n} is the total N from the same fit
#'   as h2.
#' @param gcov,intercept Stage 1-2 inputs, or \code{NULL} to estimate.
#' @param ldscore LD scores for stage 1-2 estimation.
#' @param tag_r2 Clump threshold (\code{NULL} for per-SNP estimates).
#' @param method,rectify Passed to [lrcq()].
#' @param windows Optional precomputed windows (see [ld_from_windows()]).
#' @param ... Passed to [lrcq()].
#' @return List with \code{stage12} (gcov, intercept) and \code{w}
#'   (data frame from [lrcq()], with SNP metadata when available).
#' @export
run_lrcq <- function(Z, ld, n, M, gcov = NULL, intercept = NULL, ldscore = NULL,
                     tag_r2 = 0.5, method = "wls", rectify = "C",
                     windows = NULL, ...) {
  if (is.null(gcov)) {
    if (is.null(ldscore)) stop("give gcov/intercept or ldscore for stage 1-2")
    st <- ldsc_matrix(Z, ldscore, n, M)
    gcov <- st$gcov
    if (is.null(intercept)) intercept <- st$intercept
  } else check_gcov(gcov, ncol(Z), intercept)
  if (is.null(intercept)) intercept <- diag(ncol(Z))
  w <- lrcq(Z, ld, n, gcov, M, intercept, method = method, rectify = rectify,
            tag_r2 = tag_r2, windows = windows, ...)
  if (!is.null(ld$snp)) w <- cbind(ld$snp, w)
  list(stage12 = list(gcov = gcov, intercept = intercept), w = w)
}

# Warn when a diagonal genetic covariance is supplied for several traits:
# LRCP and LRCQ standard errors then assume genetically independent traits,
# which on real phenome-wide data (e.g. many anthropometric traits) makes
# them far too small.
check_gcov <- function(gcov, q, intercept = NULL) {
  if (q > 1 && !is.null(dim(gcov)) && all(gcov[upper.tri(gcov)] == 0))
    warning("gcov is diagonal: standard errors assume genetically uncorrelated traits. ",
            "Supply genetic covariances (e.g. ldsc_matrix() or gcov_from_rg()).", call. = FALSE)
  # a common slip: univariate LDSC intercepts on the diagonal with cross-trait
  # intercepts from a joint fit off it, which can leave C indefinite
  if (q > 1 && !is.null(dim(intercept))) {
    ev <- eigen(intercept, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8 * max(ev))
      warning("intercept matrix is not positive semi-definite (smallest eigenvalue ",
              signif(min(ev), 3), "); take its diagonal and off-diagonal from the same fit.",
              call. = FALSE)
  }
  invisible(NULL)
}

#' Genetic covariance from heritabilities and genetic correlations
#'
#' @param h2 Length-q SNP heritabilities.
#' @param rg q-by-q genetic correlation matrix (projected to the nearest
#'   positive semi-definite correlation matrix if needed).
#' @return q-by-q genetic covariance matrix.
#' @export
gcov_from_rg <- function(h2, rg) {
  rg <- (rg + t(rg)) / 2
  e <- eigen(rg, symmetric = TRUE)
  if (min(e$values) < 0) {
    rg <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
    d <- sqrt(diag(rg)); rg <- rg / outer(d, d)
  }
  h <- sqrt(pmax(h2, 0))
  rg * outer(h, h)
}

#' Run LRCP for distal window pairs (stage 4)
#'
#' Screens tags on stage-3 enrichment, then estimates genetic-effect
#' correlation for every pair of screened tags in different windows that
#' share no LD (different chromosomes, or further apart than
#' \code{min_distance}). P-values are adjusted by Benjamini-Hochberg
#' (primary) and Benjamini-Yekutieli (conservative check) across all pairs.
#'
#' The default test is the conditional Gaussian test of [lrcp_gene()]
#' (theory supplement S8.2, S8.5): it uses the model for one window only and
#' is robust to noisy enrichment estimates in the other. Every pair also
#' carries the GLS-type z (\code{z_gls}), the ratio diagnostic \code{R_AB}
#' (\eqn{z/z_{gls} = R_{AB}^{-1/2}}) and the side-B signal energy ratio
#' \code{e_B}. With \code{class_null = TRUE}, screened tags are clustered
#' genome-wide by sign-free trait-profile correlation ([profile_classes()]),
#' and each pair is tested against the class null built from held-out loci
#' of each tag's class, excluding both tags' blocks ([class_covariance()]);
#' \code{q_eff_class1/2} are reported and pairs with a class effective trait
#' count below 5 are flagged \code{testable = FALSE}: report their
#' \eqn{\hat\rho} descriptively as shared-mediator alignment. Classes with
#' fewer than \code{min_class_loci} held-out loci fall back to the
#' genome-wide covariance.
#'
#' @param Z,ld,n,M As in [run_lrcq()].
#' @param gcov,intercept Stage 1-2 values. \code{gcov} must hold genome-wide
#'   LDSC genetic covariances, never a diagonal matrix.
#' @param lrcq_fit Output of [run_lrcq()] (its \code{w} table).
#' @param threshold Screening cutoff on enrichment.
#' @param z_min Optional screening cutoff on the enrichment z-score.
#' @param min_distance Minimum distance (bp) between windows on the same
#'   chromosome.
#' @param test \code{"conditional"} (default) or \code{"gls"} ([lrcp_distal()]).
#' @param method LRCP estimator for [lrcp_distal()] when \code{test = "gls"}.
#' @param denominators See [lrcp_distal()] (\code{test = "gls"}).
#' @param class_null Test against locus-class nulls (S8.5).
#' @param min_abs_cor,min_class_loci Class definition: see
#'   [profile_classes()]; minimum held-out loci per class.
#' @param n_boot Conditional bootstrap draws for the conditional test
#'   ([lrcp_gene()], S8.2a; 0 for the normal reference). Ignored with the
#'   class null, which uses the normal reference.
#' @param mle Also run [lrcp_mle()] per window pair.
#' @return Data frame with one row per tag pair; \code{p} is the bootstrap
#'   p-value (larger of both orientations) for the conditional test.
#' @export
run_lrcp <- function(Z, ld, n, M, gcov, intercept = NULL, lrcq_fit,
                     threshold = 10, z_min = -Inf, min_distance = 5e6,
                     test = c("conditional", "gls"), method = "gls",
                     denominators = "R", class_null = FALSE,
                     min_abs_cor = 0.5, min_class_loci = 5, n_boot = 1000,
                     mle = FALSE) {
  test <- match.arg(test)
  q <- ncol(Z)
  check_gcov(gcov, q, intercept)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  wt <- lrcq_fit$w
  cand <- which(!is.na(wt$w_raw) & wt$w_raw >= threshold & wt$z >= z_min)
  if (length(cand) < 2) return(data.frame())
  blk <- ld$block_id[cand]
  groups <- split(cand, blk)
  chr <- if (!is.null(ld$snp$chr)) ld$snp$chr else rep(1, ld$m)
  pos <- if (!is.null(ld$snp$pos)) ld$snp$pos else seq_len(ld$m)
  # per-block projected profiles of the screened tags
  prof <- lapply(groups, function(S) {
    Rs <- as.matrix(ld$get(S))
    v <- pmax(local_variance(Z[S, , drop = FALSE], Rs, seq_along(S), n, gcov, M, intercept), 0)
    Ri <- solve(Rs)
    list(S = S, R = Rs, v = v, X = Ri %*% Z[S, , drop = FALSE], tauR = diag(Ri))
  })
  Xall <- do.call(rbind, lapply(prof, `[[`, "X"))
  tauR <- unlist(lapply(prof, `[[`, "tauR")); tauG <- unlist(lapply(prof, `[[`, "v"))
  gid <- rep(seq_along(prof), vapply(prof, function(p) length(p$S), 1L))
  cls <- if (class_null) profile_classes(Xall, min_abs_cor) else rep(1L, nrow(Xall))
  Tm <- scale_matrix(n, gcov, M)
  vclass <- function(rows, i, j) lapply(rows, function(r) {
    use <- which(cls == cls[r] & gid != i & gid != j)
    cc <- if (length(use) >= min_class_loci) class_covariance(Xall, tauR, tauG, intercept, use) else NULL
    if (is.null(cc) || is.null(cc$V)) list(V = Tm, q_eff = sum(diag(Tm))^2 / sum(Tm^2)) else cc
  })
  out <- list()
  for (i in seq_along(groups)) for (j in seq_along(groups)) {
    if (j <= i) next
    A <- groups[[i]]; B <- groups[[j]]
    same_chr <- chr[A[1]] == chr[B[1]]
    if (same_chr && min(abs(outer(pos[A], pos[B], "-"))) < min_distance) next
    res <- data.frame(snp1 = rep(A, length(B)), snp2 = rep(B, each = length(A)))
    if (test == "conditional") {
      pa <- prof[[i]]; pb <- prof[[j]]
      VA <- VB <- NULL
      if (class_null) {
        ca <- vclass(which(gid == i), i, j); cb <- vclass(which(gid == j), i, j)
        VA <- lapply(ca, `[[`, "V"); VB <- lapply(cb, `[[`, "V")
        res$class1 <- rep(cls[gid == i], length(B)); res$class2 <- rep(cls[gid == j], each = length(A))
        res$q_eff_class1 <- rep(vapply(ca, `[[`, 0, "q_eff"), length(B))
        res$q_eff_class2 <- rep(vapply(cb, `[[`, 0, "q_eff"), each = length(A))
      }
      f <- suppressWarnings(lrcp_gene(Z[A, , drop = FALSE], Z[B, , drop = FALSE], pa$R, pb$R,
                                      n = n, gcov = gcov, M = M, intercept = intercept,
                                      VclassA = VA, VclassB = VB, n_boot = n_boot, n_sim = 0))
      den <- sqrt(outer(pa$v, pb$v))
      res$rho <- as.vector(ifelse(den > 0, f$C / den, NA_real_))
      res$C <- as.vector(f$C); res$se <- as.vector(f$se); res$z <- as.vector(f$z)
      res$p <- as.vector(f$p); res$p_maxT_window_pair <- as.vector(f$p_maxT %||% NA)
      res$z_gls <- as.vector(f$z_gls); res$R_AB <- as.vector(f$R_AB); res$e_B <- as.vector(f$e_B)
      if (class_null) res$testable <- pmin(res$q_eff_class1, res$q_eff_class2) >= 5
    } else {
      ia <- which(ld$block_id == ld$block_id[A[1]])
      ib <- which(ld$block_id == ld$block_id[B[1]])
      RA <- ld$get(ia); RB <- ld$get(ib)
      wa <- numeric(length(ia)); wa[match(A, ia)] <- wt$w_raw[A]
      wb <- numeric(length(ib)); wb[match(B, ib)] <- wt$w_raw[B]
      f <- lrcp_distal(Z[ia, , drop = FALSE], Z[ib, , drop = FALSE], RA, RB, wa, wb,
                       n, gcov, M, intercept, S1 = match(A, ia), S2 = match(B, ib),
                       method = method, denominators = denominators)
      res$rho <- as.vector(f$rho); res$se <- as.vector(f$se)
      res$z <- res$rho / res$se
    }
    if (mle) {
      ia <- which(ld$block_id == ld$block_id[A[1]])
      ib <- which(ld$block_id == ld$block_id[B[1]])
      RA <- ld$get(ia); RB <- ld$get(ib)
      wa <- numeric(length(ia)); wa[match(A, ia)] <- wt$w_raw[A]
      wb <- numeric(length(ib)); wb[match(B, ib)] <- wt$w_raw[B]
      g <- lrcp_mle(Z[ia, , drop = FALSE], Z[ib, , drop = FALSE], RA, RB, wa, wb,
                    n, gcov, M, intercept, S1 = match(A, ia), S2 = match(B, ib))
      res$rho_mle <- as.vector(g$rho); res$se_mle <- as.vector(g$se)
      res$lrt_window_pair <- g$lrt
    }
    out[[length(out) + 1]] <- res
  }
  res <- do.call(rbind, out)
  if (is.null(res)) return(data.frame())
  if (class_null && any(!res$testable))
    warning(sum(!res$testable), " pairs have a class effective trait count below 5: ",
            "report their rho descriptively, not as tests (theory S8.5)", call. = FALSE)
  if (is.null(res$p)) res$p <- 2 * stats::pnorm(-abs(res$z))
  res$q_BH <- stats::p.adjust(res$p, "BH")
  res$q_BY <- stats::p.adjust(res$p, "BY")
  if (!is.null(ld$snp)) {
    id <- if (!is.null(ld$snp$ID)) ld$snp$ID else paste(chr, pos, sep = ":")
    res$id1 <- id[res$snp1]; res$id2 <- id[res$snp2]
  }
  res
}

#' Choose genetically diverse traits
#'
#' Greedy pruning on genetic correlation: traits are taken in order of
#' \code{priority} (e.g. h2 z-score) and kept only if their absolute genetic
#' correlation with every trait already kept is at most \code{max_rg}. Power
#' of LRCQ and LRCP grows with the number of genetically independent trait
#' dimensions, not with the number of traits (theory supplement S9).
#'
#' @param rg q-by-q genetic correlation matrix.
#' @param priority Length-q score; higher is taken first.
#' @param max_rg Largest absolute genetic correlation allowed between kept
#'   traits.
#' @param max_traits Optional cap on the number of traits kept.
#' @return Indices of kept traits, in the order chosen.
#' @export
prune_traits <- function(rg, priority = rep(1, nrow(rg)), max_rg = 0.5,
                         max_traits = Inf) {
  ord <- order(-priority)
  keep <- integer(0)
  for (i in ord) {
    if (length(keep) >= max_traits) break
    if (!length(keep) || all(abs(rg[i, keep]) <= max_rg | is.na(rg[i, keep]))) keep <- c(keep, i)
  }
  keep
}

#' Check the scale of n, h2 and M against the data
#'
#' Per-trait ratio of the observed mean chi-square signal to the signal
#' expected under \eqn{w = 1}: \eqn{(\overline{z^2} - C_{aa}) /
#' (n_a h^2_a \bar\ell / M)}, with \eqn{\ell} the LD scores of the analysed
#' SNPs computed within the same reference as \code{R}. A median near 1 means
#' n, h2 and M are on one scale. \code{M} is the number of SNPs the model
#' spreads h2 over, i.e. the SNPs whose LD enters \code{R} (1,094,844 for the
#' Pan-UKB HapMap3 set), not the M_5_50 count that Pan-UKB's own LDSC used with
#' full-panel LD scores. In the TP53/MDM2 smoke test, Pan-UKB h2 with total N gave
#' a median of 0.95 with M = 1,094,844 and 5.9 with M_5_50 = 6,805,960. Use
#' genome-wide or large regions; small regions are noisy and enriched.
#'
#' @param Z,ld,n,M As in [lrcq()].
#' @param h2 Length-q heritabilities from the same fit as \code{n}.
#' @param intercept Length-q univariate LDSC intercepts (default 1).
#' @return Length-q vector of ratios.
#' @export
check_scale <- function(Z, ld, n, h2, M, intercept = 1) {
  Z <- as.matrix(Z)
  l <- mean(ld_scores(ld))
  intercept <- rep_len(intercept, ncol(Z))
  (colMeans(Z^2) - intercept) / (n * h2 * l / M)
}
