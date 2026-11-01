test_that("command-line munge and lrcq run on the bundled example", {
  ex <- system.file("extdata", "example", package = "lrcpq")
  skip_if(ex == "")
  td <- tempfile(); dir.create(td)
  ids <- sprintf("trait%02d", 1:3)
  for (i in ids)
    lrcpq_cli(c("munge", "--sumstats", file.path(ex, paste0(i, ".sumstats.gz")),
                "--ref-ld", file.path(ex, "ref_chr@"), "--chr", "22",
                "--out", file.path(td, paste0(i, ".munged.gz"))))
  m <- utils::read.delim(file.path(td, "trait01.munged.gz"))
  expect_true(all(c("chr", "pos", "rsid", "a1", "a2", "z", "n") %in% names(m)))
  expect_gt(nrow(m), 1900)
  writeLines(paste0(ids, "\t", file.path(td, paste0(ids, ".munged.gz"))), file.path(td, "list.txt"))
  h2 <- utils::read.table(file.path(ex, "h2.tsv"))
  rg <- as.matrix(utils::read.delim(file.path(ex, "rg.tsv"), row.names = 1))[ids, ids]
  utils::write.table(h2[h2$V1 %in% ids, ], file.path(td, "h2.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(rg, file.path(td, "rg.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  suppressWarnings(lrcpq_cli(c("lrcq", "--sumstats-list", file.path(td, "list.txt"),
              "--ref-ld", file.path(ex, "ref_chr@"), "--chr", "22", "--M", "1e6",
              "--h2", file.path(td, "h2.tsv"), "--rg", file.path(td, "rg.tsv"),
              "--method", "ols", "--out", file.path(td, "ex"))))
  w <- utils::read.delim(file.path(td, "ex.lrcq.tsv.gz"))
  expect_true(all(c("rsid", "w_raw", "se", "w_null1") %in% names(w)))
  expect_true(file.exists(file.path(td, "ex.log")))
})
