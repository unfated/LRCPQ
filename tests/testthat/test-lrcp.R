sim_two_regions <- function(q = 30, rho = c(0.6, -0.4, 0), gcov = NULL, C = NULL) {
  blocks <- make_ld_ar1(c(30, 30), 0.5); m <- 60
  w <- numeric(m); idx1 <- c(5, 12, 20); idx2 <- 30 + c(8, 15, 25); w[c(idx1, idx2)] <- 10
  Rb <- as.matrix(make_Rb(m, "pairs", pairs = cbind(idx1, idx2), pair_rho = rho))
  if (is.null(gcov)) gcov <- make_gcov(q, h2 = rep(0.3, q))$gcov
  if (is.null(C)) C <- diag(q)
  list(blocks = blocks, w = w, Rb = Rb, gcov = gcov, C = C, n = rep(3e5, q), M = 2e4)
}

test_that("distal LRCP regression is unbiased, calibrated, and halves under the chapter form", {
  set.seed(11)
  s <- sim_two_regions(q = 30)
  est <- se <- chap <- NULL
  for (r in 1:120) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, Rb = s$Rb, gcov = s$gcov, M = s$M)
    a <- list(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]], s$w[1:30], s$w[31:60], s$n, s$gcov, s$M)
    f <- do.call(lrcp_distal, a)
    g <- do.call(lrcp_distal, c(a, list(pair_sum = "chapter")))
    est <- rbind(est, diag(f$rho)); se <- rbind(se, diag(f$se)); chap <- rbind(chap, diag(g$rho))
  }
  sdv <- apply(est, 2, sd)
  expect_true(all(abs(colMeans(est) - c(0.6, -0.4, 0)) < 4 * sdv / sqrt(120)))
  expect_equal(colMeans(chap), colMeans(est) / 2, tolerance = 1e-8)
  expect_true(all(abs(colMeans(se) / sdv - 1) < 0.25))
})

test_that("LRCP SEs account for correlated traits and sample overlap", {
  set.seed(12)
  q <- 30
  gc <- make_gcov(q, h2 = rep(0.3, q), type = "cluster", rg = 0.7, n_clusters = 2)
  s <- sim_two_regions(q, rho = c(0, 0, 0), gcov = gc$gcov, C = make_intercept(q, 1, 0.4))
  zz <- NULL
  for (r in 1:120) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, gcov = s$gcov, intercept = s$C, M = s$M)
    f <- lrcp_distal(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                     s$w[1:30], s$w[31:60], s$n, s$gcov, s$M, s$C)
    zz <- c(zz, as.vector(f$z))
  }
  expect_lt(abs(sd(zz) - 1), 0.15)
})

test_that("LRCP MLE recovers rho and the gradient is exact", {
  set.seed(13)
  s <- sim_two_regions(q = 30)
  est <- NULL
  for (r in 1:15) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, Rb = s$Rb, gcov = s$gcov, M = s$M)
    f <- lrcp_mle(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                  s$w[1:30], s$w[31:60], s$n, s$gcov, s$M)
    expect_equal(f$convergence, 0)
    expect_gt(f$lrt, -1e-6)
    est <- rbind(est, diag(f$rho))
  }
  expect_true(all(abs(colMeans(est) - c(0.6, -0.4, 0)) < 0.12))
})

test_that("trait covariance helper matches the moment equation", {
  R <- make_ld_ar1(5, 0.5)[[1]]; w <- c(2, 0, 1, 1, 1)
  g <- make_gcov(3, h2 = c(0.2, 0.3, 0.4), type = "cs", rg = 0.5)$gcov
  n <- c(1e4, 2e4, 3e4); C <- make_intercept(3, 0.5, 0.2)
  Mk <- expected_trait_cov(2, R, w, NULL, n, g, 100, C)
  for (a in 1:3) for (b in 1:3)
    expect_equal(Mk[a, b], expected_zz(R, w, NULL, n, g, 100, C, a, b)[2, 2])
  expect_equal(q_eff(diag(3)), 3)
})

test_that("R-projected local variance estimates w of candidates", {
  set.seed(23)
  s <- sim_two_regions(q = 60)
  v <- NULL
  for (r in 1:40) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, gcov = s$gcov, M = s$M)
    v <- rbind(v, local_variance(sim$Z[1:30, ], s$blocks[[1]], c(5, 12, 20), s$n, s$gcov, s$M))
  }
  expect_true(all(abs(colMeans(v) - 10) < 4 * apply(v, 2, sd) / sqrt(40) + 0.5))
})

test_that("pipelines run end to end and rank true pairs first", {
  set.seed(61)
  blocks <- make_ld_ar1(rep(50, 6), 0.7); ld <- ld_from_blocks(blocks); m <- ld$m
  w <- numeric(m); cand <- c(10, 60, 110, 160, 210, 260); w[cand] <- 40
  w[sample(setdiff(1:m, cand), 30)] <- 2; w <- w * m / sum(w)
  Rb <- as.matrix(make_Rb(m, "pairs", pairs = rbind(c(10, 160), c(60, 210)), pair_rho = c(0.7, -0.6)))
  q <- 50; gc <- make_gcov(q, h2 = rep(0.3, q), type = "cluster", rg = 0.4, n_clusters = 5)
  C <- make_intercept(q, 1, 0.2); n <- rep(4e5, q); M <- 1e5
  sim <- simulate_lrcpq(blocks, n = n, w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = M)
  fq <- run_lrcq(sim$Z, ld, n, M, gc$gcov, C)
  expect_true(all(cand %in% which(fq$w$w_raw > 5)))
  fp <- run_lrcp(sim$Z, ld, n, M, gc$gcov, C, fq, threshold = 5, min_distance = 0)
  top <- fp[order(fp$p), ][1:2, ]
  expect_setequal(paste(top$snp1, top$snp2), c("10 160", "60 210"))
})

test_that("plug-in LRCP SE is calibrated at large rho", {
  set.seed(14)
  s <- sim_two_regions(q = 30, rho = c(0.9, 0.3, 0))
  est <- se <- NULL
  for (r in 1:150) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, Rb = s$Rb, gcov = s$gcov, M = s$M)
    f <- lrcp_distal(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                     s$w[1:30], s$w[31:60], s$n, s$gcov, s$M, se_type = "plugin")
    est <- rbind(est, diag(f$rho)); se <- rbind(se, diag(f$se))
  }
  expect_lt(abs(mean(se[, 1]) / sd(est[, 1]) - 1), 0.15)
})

test_that("gcov_from_rg builds a PSD covariance and diagonal gcov warns", {
  rg <- matrix(c(1, 0.9, -0.9, 0.9, 1, 0.9, -0.9, 0.9, 1), 3)   # not PSD
  G <- gcov_from_rg(c(0.2, 0.3, 0.4), rg)
  expect_gte(min(eigen(G)$values), -1e-10)
  expect_equal(diag(G), c(0.2, 0.3, 0.4))
  expect_warning(lrcpq:::check_gcov(diag(3), 3), "diagonal")
  C <- matrix(c(1.2, 1.4, 1.4, 1.2), 2)
  expect_warning(lrcpq:::check_gcov(matrix(0.1, 2, 2), 2, C), "positive semi-definite")
  expect_silent(lrcpq:::check_gcov(matrix(0.1, 2, 2), 2, diag(2)))
})

test_that("prune_traits keeps a genetically diverse set", {
  rg <- matrix(0.1, 4, 4); diag(rg) <- 1; rg[1, 2] <- rg[2, 1] <- 0.9
  expect_equal(prune_traits(rg, c(1, 2, 0.5, 0.2)), c(2, 3, 4))
  expect_equal(prune_traits(rg, c(1, 2, 0.5, 0.2), max_traits = 2), c(2, 3))
})

test_that("candidates with non-positive variance are dropped, not an error", {
  set.seed(15)
  s <- sim_two_regions(q = 30)
  sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, gcov = s$gcov, M = s$M)
  w1 <- s$w[1:30]; w1[c(5, 2)] <- c(10, 0)
  for (m in c("gls", "ols", "wls")) {
    f <- suppressWarnings(lrcp_distal(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                                      w1, s$w[31:60], s$n, s$gcov, s$M, S1 = c(2, 5, 12, 20),
                                      S2 = which(s$w[31:60] > 0), method = m))
    expect_true(all(is.na(f$rho[1, ])))
    expect_true(all(is.finite(f$rho[-1, ])))
  }
  expect_warning(lrcp_distal(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                             w1, s$w[31:60], s$n, s$gcov, s$M, S1 = c(2, 5, 12, 20)), "dropped")
})

test_that("check_scale is about 1 when n, h2 and M match the simulation", {
  set.seed(5)
  b <- make_ld_ar1(rep(50, 20), 0.5); ld <- ld_from_blocks(b)
  q <- 5; n <- rep(2e5, q); g <- make_gcov(q, h2 = rep(0.3, q))$gcov
  sim <- simulate_lrcpq(b, n = n, w = rep(1, ld$m), gcov = g, M = 2e4)
  r <- check_scale(sim$Z, ld, n, diag(g), M = 2e4)
  expect_true(abs(median(r) - 1) < 0.25)
  expect_true(median(check_scale(sim$Z, ld, n, diag(g), M = 2e5)) > 5)
})

test_that("lrcp_screen: simulated-null calibration is valid under H0", {
  set.seed(8)
  q <- 30; n <- rep(3e5, q); M <- 2e4
  gc <- make_gcov(q, type = "cluster", rg = 0.4); C <- make_intercept(q, overlap = 1, rp = 0.2)
  b <- make_ld_ar1(rep(4, 12), 0.4)
  w <- rep(c(30, 10, 5, 1), 12)
  sim <- simulate_lrcpq(b, n = n, w = w, gcov = gc$gcov, intercept = C, M = M)
  win <- rep(seq_along(b), each = 4); names(b) <- seq_along(b)
  sc <- lrcp_screen(sim$Z, win, b, chr = win, pos = rep(1, 48), n = n, gcov = gc$gcov, M = M,
                    intercept = C, n_null = 20, p_report = 1)
  expect_equal(sc$summary$n_pairs, choose(48, 2) - 12 * choose(4, 2))
  expect_lt(sc$summary$null_sd, 1)
  expect_lt(abs(sc$summary$sd_z - sc$summary$null_sd), 0.15)
  expect_lt(mean(sc$pairs$p_cal < 0.05), 0.09)
  expect_true(all(sc$pairs$p_cal >= sc$pairs$p_norm - 0.2))
})
