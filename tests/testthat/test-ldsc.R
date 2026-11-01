test_that("LDSC recovers heritability, co-heritability and intercepts", {
  set.seed(5)
  blocks <- make_ld_ar1(rep(50, 40), 0.7)
  ld <- ld_from_blocks(blocks)
  m <- ld$m; q <- 2
  g <- make_gcov(q, h2 = c(0.5, 0.3), type = "cs", rg = 0.5)$gcov
  C <- make_intercept(q, overlap = 1, rp = 0.3, diag_intercept = 1)
  l <- ld_scores(ld)
  res <- replicate(30, {
    sim <- simulate_lrcpq(blocks, n = c(5e3, 5e3), gcov = g, intercept = C, M = m, prop_nonzero = 1, dist = "point")
    f <- ldsc_matrix(sim$Z, l, sim$n, m, n_blocks = 50)
    c(f$h2, f$gcov[1, 2], f$intercept[1, 1], f$intercept[1, 2])
  })
  est <- rowMeans(res); sd <- apply(res, 1, sd) / sqrt(30)
  truth <- c(0.5, 0.3, g[1, 2], 1, 0.3)
  expect_true(all(abs(est - truth) < 4 * sd + 0.01))
})
