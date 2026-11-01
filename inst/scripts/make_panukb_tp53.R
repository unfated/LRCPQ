# Builds inst/extdata/panukb_tp53 for the real-data vignette: Pan-UKB EUR
# Z-scores (CC BY 4.0) of 20 traits at 2,007 HapMap3 SNPs around TP53
# (chr17, three EUR LD blocks), aligned to A1 of the 1000G EUR HM3 reference,
# with the genome-wide LDSC rg and intercepts and Pan-UKB h2 of the 444-trait
# analysis. Inputs: the fetch (panukb_z) saved by the software thread, the
# 444-trait LDSC matrices, and the 1000G reference (GitHub release asset).
# Run from the package root: Rscript inst/scripts/make_panukb_tp53.R FETCH.rds REALDIR REFDIR PLINK
library(lrcpq)
a <- commandArgs(TRUE); f <- readRDS(a[1]); real <- a[2]; refdir <- a[3]; plink <- a[4]
out <- "inst/extdata/panukb_tp53"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
r <- f$ref; k <- nrow(r); Zall <- f$Zall
Z <- ifelse(!is.na(Zall[1:k, ]), Zall[1:k, ], -Zall[k + 1:k, ])      # aligned to A1
ok <- rowSums(is.na(Z)) == 0 & !paste0(r$a1, r$a2) %in% c("AT", "TA", "CG", "GC")
ids <- f$pick$trait_id
colnames(Z) <- ids
h2t <- read.delim(file.path(real, "panukb_ldsc.h2.tsv"))
traits <- data.frame(trait_id = ids, description = f$pick$description, N = f$pick$N,
                     h2 = h2t$h2_panukb_ldsc[match(ids, h2t$trait_id)])
write.table(traits, file.path(out, "traits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
d <- data.frame(chr = r$chr, pos = r$pos, rsid = r$rsid, a1 = r$a1, a2 = r$a2, signif(Z, 5),
                check.names = FALSE)[ok, ]
con <- gzfile(file.path(out, "z.tsv.gz"), "w"); write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
for (m in c("rg", "intercept")) {
  x <- as.matrix(read.delim(file.path(real, paste0("panukb_ldsc.", m, ".tsv")), row.names = 1, check.names = FALSE))
  write.table(signif(x[ids, ids], 6), file.path(out, paste0(m, ".tsv")), sep = "\t", quote = FALSE, col.names = NA)
}
system2(plink, c("--bfile", file.path(refdir, "eur_hm3_chr17"), "--chr", "17", "--from-bp", f$region$start,
                 "--to-bp", f$region$end - 1, "--make-bed", "--out", file.path(out, "ref_chr17_tp53")))
unlink(file.path(out, c("ref_chr17_tp53.log", "ref_chr17_tp53.nosex")))
cat(sum(ok), "SNPs x", length(ids), "traits written to", out, "\n")
