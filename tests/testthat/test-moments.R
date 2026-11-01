test_that("matrix and element forms of the moment equation agree", {
  set.seed(1)
  m <- 5; q <- 2
  R <- make_ld_ar1(m, 0.6)[[1]]
  w <- make_w(m, 0.6)
  Rb <- as.matrix(make_Rb(m, "pairs", pairs = rbind(c(1, 4), c(2, 5)), pair_rho = c(0.5, -0.3)))
  g <- make_gcov(q, h2 = c(0.3, 0.2), type = "cs", rg = 0.4)$gcov
  n <- c(1e4, 2e4); C <- make_intercept(q, 0.5, 0.2)
  E_once <- expected_zz(R, w, Rb, n, g, M = m, intercept = C, a = 1, b = 2)
  E_chap <- expected_zz(R, w, Rb, n, g, M = m, intercept = C, a = 1, b = 2, pair_sum = "chapter")
  for (k in 1:m) for (l in 1:m) {
    expect_equal(E_once[k, l], expected_cross_element(k, l, 1, 2, R, w, Rb, n, g, m, C, 0.5))
    expect_equal(E_chap[k, l], expected_cross_element(k, l, 1, 2, R, w, Rb, n, g, m, C, 1))
  }
  E2 <- expected_z2(R, w, Rb, n, g, M = m, intercept = C)
  expect_equal(E2[, 1], diag(expected_zz(R, w, Rb, n, g, m, C, 1, 1)))
})

test_that("simulated Z-scores match the expected moments", {
  set.seed(2)
  m <- 6; q <- 3
  R <- make_ld_ar1(m, 0.5)[[1]]
  w <- make_w(m, 0.5)
  Rb <- as.matrix(make_Rb(m, "cluster", clusters = list(which(w > 0)), rho = 0.6))
  g <- make_gcov(q, h2 = c(0.5, 0.4, 0.3), type = "cs", rg = 0.5)$gcov
  n <- c(100, 200, 150); C <- make_intercept(q, 1, 0.3)
  reps <- 20000
  acc12 <- matrix(0, m, m); acc11 <- matrix(0, m, m)
  for (r in 1:reps) {
    B <- simulate_B(w, Rb, g, M = m)
    Z <- simulate_Z(R, B, n, C)
    acc12 <- acc12 + tcrossprod(Z[, 1], Z[, 2])
    acc11 <- acc11 + tcrossprod(Z[, 1])
  }
  E12 <- expected_zz(R, w, Rb, n, g, m, C, 1, 2)
  E11 <- expected_zz(R, w, Rb, n, g, m, C, 1, 1)
  expect_lt(max(abs(acc12 / reps - E12) / (abs(E11) + 1)), 0.06)
  expect_lt(max(abs(acc11 / reps - E11) / (abs(E11) + 1)), 0.06)
})
