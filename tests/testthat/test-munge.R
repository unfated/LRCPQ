test_that("munge_sumstats aligns to A1, converts BETA/SE and filters", {
  ref <- data.frame(chr = c(1, 1, 1, 2), pos = c(10, 20, 30, 40), a1 = c("A", "C", "A", "G"),
                    a2 = c("G", "T", "T", "A"), rsid = paste0("rs", 1:4))
  d <- data.frame(CHR = c("chr1", "1", "1", "2", "2"), BP = c(10, 20, 30, 40, 50),
                  ALT = c("A", "T", "A", "G", "C"), REF = c("G", "C", "T", "A", "G"),
                  BETA = c(0.2, 0.3, 0.1, -0.4, 1), SE = 0.1, N_CASES = 1000, N_CONTROLS = 3000,
                  INFO = c(1, 1, 1, 0.5, 1))
  m <- munge_sumstats(d, ref)
  expect_equal(m$j, 1:2)                     # rs3 ambiguous, rs4 low INFO, pos 50 absent
  expect_equal(m$z, c(2, -3))                # second SNP flipped to the reference A1
  expect_equal(m$n, c(4000, 4000))           # total N
  zz <- build_z_matrix(list(t1 = m, t2 = m), nrow(ref))
  expect_equal(zz$keep, 1:2); expect_equal(unname(zz$n), c(4000, 4000))
  expect_error(munge_sumstats(d[, c("CHR", "BP", "BETA", "SE")], ref), "allele")
})
