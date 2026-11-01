test_that("fast OLS/WLS equal the brute-force stacked regression", {
  set.seed(3)
  m <- 8; q <- 4
  R <- make_ld_ar1(m, 0.5)[[1]]
  sim <- simulate_lrcpq(R, n = c(5e3, 1e4, 2e4, 8e3), q = q, M = m, prop_nonzero = 0.5)
  f_ols <- lrcq_window(sim$Z, R, sim$n, sim$gcov, sim$M, method = "ols", se = "none")
  b_ols <- lrcq_bruteforce(sim$Z, R, sim$n, sim$gcov, sim$M)
  expect_equal(f_ols$w, unname(b_ols), tolerance = 1e-8)
  # WLS with weights from the OLS fit
  s <- sim$n * diag(sim$gcov) / m
  V <- 2 * (outer(pmax(as.vector(R^2 %*% f_ols$w), 0), s) + 1)^2
  f_wls <- lrcq_window(sim$Z, R, sim$n, sim$gcov, sim$M, method = "wls", se = "none")
  b_wls <- lrcq_bruteforce(sim$Z, R, sim$n, sim$gcov, sim$M, V = V)
  expect_equal(f_wls$w, unname(b_wls), tolerance = 1e-8)
  # equal weighting is the average of per-trait D^{-1} y_a / s_a
  f_eq <- lrcq_window(sim$Z, R, sim$n, sim$gcov, sim$M, method = "equal", se = "none")
  Y <- sim$Z^2 - 1
  manual <- rowMeans(solve(R^2, Y) / matrix(s, m, q, byrow = TRUE))
  expect_equal(f_eq$w, manual, tolerance = 1e-8)
})

test_that("LRCQ is unbiased with known parameters and SEs are calibrated", {
  set.seed(4)
  blocks <- make_ld_ar1(rep(20, 3), 0.4)
  ld <- ld_from_blocks(blocks)
  m <- ld$m; q <- 20
  w <- make_w(m, 0.3)
  g <- make_gcov(q, h2 = rep(0.4, q))$gcov
  n <- rep(2e4, q)
  est <- ses <- NULL
  for (r in 1:150) {
    sim <- simulate_lrcpq(blocks, n = n, w = w, gcov = g, M = m)
    f <- lrcq(sim$Z, ld, n, g, M = m, method = "wls")
    est <- rbind(est, f$w_raw); ses <- rbind(ses, f$se)
  }
  bias <- colMeans(est) - w
  expect_lt(max(abs(bias) / (apply(est, 2, sd) / sqrt(150))), 4.5)
  ratio <- colMeans(ses) / apply(est, 2, sd)
  expect_gt(mean(ratio), 0.85); expect_lt(mean(ratio), 1.15)
})

test_that("rectification methods follow chapter 6.3.2", {
  w <- c(-0.5, -0.2, 0.1, 0.3, 0.6, 2, 5)
  expect_equal(rectify_w(w, "A"), c(0, 0, 0, 0, 0.6, 2, 5))
  expect_equal(rectify_w(w, "B"), c(0, 0, 0, 0, 0.6, 2, 5))
  expect_equal(rectify_w(w, "C"), c(0, 0, 0, 0, 0.6, 2, 5))
  w2 <- c(-0.1, 0.05, 0.3, 1)
  expect_equal(rectify_w(w2, "A"), c(0, 0, 0.3, 1))
  expect_equal(rectify_w(w2, "B"), c(0, 0, 0.3, 1))
  expect_equal(rectify_w(w2, "C"), c(0, 0, 0.3, 1))
  expect_equal(rectify_w(c(1, 2), "A"), c(1, 2))
})

test_that("tag-set LRCQ returns one calibrated estimate per tag in tight LD", {
  set.seed(22)
  blocks <- make_ld_ar1(rep(40, 3), 0.95); ld <- ld_from_blocks(blocks); m <- ld$m
  w <- make_w(m, 0.2); q <- 20; g <- make_gcov(q, h2 = rep(0.3, q))$gcov; n <- rep(3e5, q)
  est <- ses <- NULL
  for (r in 1:60) {
    sim <- simulate_lrcpq(blocks, n = n, w = w, gcov = g, M = 2e4)
    f <- lrcq(sim$Z, ld, n, g, 2e4, method = "ols", tag_r2 = 0.5)
    est <- rbind(est, f$w_raw); ses <- rbind(ses, f$se)
  }
  tg <- which(!is.na(est[1, ]))
  expect_lt(length(tg), m / 3)
  expect_true(all(f$tag[tg] == tg))
  ratio <- colMeans(ses[, tg]) / apply(est[, tg], 2, sd)
  expect_gt(median(ratio), 0.8); expect_lt(median(ratio), 1.25)
})

test_that("skew-corrected LRCQ: third cumulant matches Monte Carlo", {
  set.seed(4)
  R <- make_ld_ar1(30, 0.6)[[1]]; q <- 10
  gc <- make_gcov(q, type = "cluster", rg = 0.6); C <- make_intercept(q, overlap = 1, rp = 0.3)
  w <- make_w(30, prop_nonzero = 0.3, dist = "lnorm"); n <- rep(3e5, q); M <- 2e4
  sim <- simulate_lrcpq(list(R), n = n, w = w, gcov = gc$gcov, intercept = C, M = M)
  P <- prune_ld(R, 0.5)
  f <- lrcq_window(sim$Z, R, n, gc$gcov, M, C, method = "ols", tags = P, skew = TRUE,
                   skew_tol = 0, w_plugin = w[P])
  ws <- replicate(1500, lrcq_window(simulate_lrcpq(list(R), n = n, w = w, gcov = gc$gcov,
    intercept = C, M = M)$Z, R, n, gc$gcov, M, C, method = "ols", tags = P, se = "none")$w)
  mc <- apply(ws, 1, function(x) mean((x - mean(x))^3))
  expect_gt(cor(f$k3, mc), 0.95)
  expect_error(lrcq_window(sim$Z, R, n, gc$gcov, M, C, method = "wls", skew = TRUE), "ols")
  expect_equal(p_skew(1.96, 0, 1, 0), stats::pnorm(-1.96))
  expect_gt(p_skew(4, 0, 1, 1), stats::pnorm(-4))   # right skew: heavier upper tail
})

test_that("default SE plug-in is background-filled rule A (LRCQ paper, S3.10a)", {
  set.seed(4); R <- make_ld_ar1(40, 0.8)[[1]]; q <- 10
  gc <- make_gcov(q, type = "cluster", rg = 0.4); C <- make_intercept(q, overlap = 1, rp = 0.2)
  n <- rep(3e5, q); M <- 2e4
  sim <- simulate_lrcpq(R, n = n, w = make_w(40, 0.5), gcov = gc$gcov, intercept = C, M = M)
  P <- prune_ld(R, 0.5)
  a <- lrcq_window(sim$Z, R, n, gc$gcov, M, C, method = "ols", tags = P)
  fill_bg <- function(w) { keep <- w > abs(min(w)); bg <- if (any(!keep)) max(mean(w[!keep]), 0) else 0; ifelse(keep, w, bg) }
  b <- lrcq_window(sim$Z, R, n, gc$gcov, M, C, method = "ols", tags = P, w_plugin = fill_bg(a$w))
  expect_equal(a$se, b$se, tolerance = 1e-12)
  old <- lrcq_window(sim$Z, R, n, gc$gcov, M, C, method = "ols", tags = P, se_rectify = "A")
  expect_true(all(a$se >= old$se - 1e-12))
})
