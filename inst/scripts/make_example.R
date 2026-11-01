# Builds the bundled example (inst/extdata/example): 50 simulated GWAS on two
# 2.5 Mb regions of chromosome 22 (14 Mb apart) with real 1000G EUR LD (503
# samples, HapMap3 SNPs). Truth: 30 enriched SNPs (w = 40) and two distal
# pairs with genetic-effect correlation 0.9 and -0.9. Files are written in a
# raw GWAS format (CHR BP SNP A1 A2 BETA SE N, A1 = effect allele, about half
# of the alleles swapped) so that the tutorial starts from munging.
# Run from the package root: Rscript inst/scripts/make_example.R
library(lrcpq); set.seed(2026)
dir <- "inst/extdata/example"
ld <- ld_from_plink(c("22" = file.path(dir, "ref_chr22")))
blocks <- lapply(split(seq_len(ld$m), ld$block_id), ld$get)
m <- ld$m; q <- 50; M <- 1e6
reg2 <- ld$snp$pos >= 3e7
w <- rep(1, m)
en1 <- sample(which(!reg2), 15); en2 <- sample(which(reg2), 15)
w[c(en1, en2)] <- 300
Rb <- make_Rb(m, "pairs", pairs = cbind(en1[1:2], en2[1:2]), pair_rho = c(0.9, -0.9))
gc <- make_gcov(q, type = "cluster", rg = 0.2, n_clusters = 5)
C <- make_intercept(q, overlap = 1, rp = 0.15)
n <- round(runif(q, 3e5, 5e5))
sim <- simulate_lrcpq(blocks, n = n, w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = M)
ids <- sprintf("trait%02d", seq_len(q))
for (a in seq_len(q)) {
  maf <- ld$snp$maf; se <- 1 / sqrt(2 * maf * (1 - maf) * n[a])
  flip <- runif(m) < 0.5
  d <- data.frame(CHR = 22, BP = ld$snp$pos, SNP = ld$snp$rsid,
                  A1 = ifelse(flip, ld$snp$a2, ld$snp$a1), A2 = ifelse(flip, ld$snp$a1, ld$snp$a2),
                  BETA = signif(ifelse(flip, -1, 1) * sim$Z[, a] * se, 5), SE = signif(se, 5), N = n[a])
  con <- gzfile(file.path(dir, paste0(ids[a], ".sumstats.gz")), "w")
  write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
}
wm <- function(x, f) { dimnames(x) <- list(ids, ids)
  write.table(signif(x, 6), file.path(dir, f), sep = "\t", quote = FALSE, col.names = NA) }
h2 <- diag(gc$gcov); rg <- cov2cor(gc$gcov)
write.table(data.frame(ids, signif(h2, 6)), file.path(dir, "h2.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
wm(rg, "rg.tsv"); wm(C, "intercept.tsv")
truth <- data.frame(rsid = ld$snp$rsid, w = w)
write.table(truth[w > 1, ], file.path(dir, "truth_enriched.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(rsid1 = ld$snp$rsid[en1[1:2]], rsid2 = ld$snp$rsid[en2[1:2]], rho = c(0.9, -0.9)),
            file.path(dir, "truth_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("M used for the simulation:", M, "\n")
