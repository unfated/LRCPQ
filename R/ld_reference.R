# Real LD references: PLINK reader, LD blocks, allele harmonisation.

#' Berisa and Pickrell (2016) EUR LD blocks (GRCh37)
#'
#' Approximately independent LD blocks from LDetect, shipped with the
#' package (1,703 blocks).
#' @return Data frame with \code{chr} (integer), \code{start}, \code{stop}.
#' @export
ld_blocks_eur <- function() {
  f <- system.file("extdata", "ldetect_EUR_hg19.bed", package = "lrcpq")
  b <- utils::read.table(f, header = TRUE, strip.white = TRUE,
                         stringsAsFactors = FALSE)
  data.frame(chr = as.integer(sub("chr", "", b$chr)), start = b$start,
             stop = b$stop)
}

#' Read a PLINK .bim file
#' @param prefix PLINK file prefix.
#' @export
read_bim <- function(prefix) {
  b <- utils::read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE,
                         colClasses = c("integer", "character", "numeric",
                                        "integer", "character", "character"))
  names(b) <- c("chr", "rsid", "cm", "pos", "a1", "a2")
  b
}

#' Read genotypes from a PLINK .bed file
#'
#' Pure-R reader for SNP-major PLINK 1 binary files.
#' @param prefix PLINK file prefix.
#' @param snps Integer indices of SNPs (rows of the .bim) to read; default all.
#' @param n Number of individuals (read from the .fam when \code{NULL}).
#' @return Individuals-by-SNPs matrix of A1 allele counts (\code{NA} missing).
#' @export
read_bed <- function(prefix, snps = NULL, n = NULL) {
  if (is.null(n)) n <- length(readLines(paste0(prefix, ".fam")))
  nbytes <- ceiling(n / 4)
  if (is.null(snps)) {
    nsnp <- length(readLines(paste0(prefix, ".bim")))
    snps <- seq_len(nsnp)
  }
  con <- file(paste0(prefix, ".bed"), "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  if (!identical(magic, as.raw(c(0x6c, 0x1b, 0x01))))
    stop("not a SNP-major PLINK .bed file")
  # read the contiguous span once, then decode
  lo <- min(snps); hi <- max(snps)
  seek(con, 3 + (lo - 1) * nbytes)
  raw <- readBin(con, "raw", (hi - lo + 1) * nbytes)
  raw <- matrix(raw, nbytes, hi - lo + 1)[, snps - lo + 1, drop = FALSE]
  v <- as.integer(raw)
  codes <- rbind(bitwAnd(v, 3L), bitwAnd(bitwShiftR(v, 2L), 3L),
                 bitwAnd(bitwShiftR(v, 4L), 3L), bitwShiftR(v, 6L))
  # 00 hom A1 -> 2, 01 missing, 10 het -> 1, 11 hom A2 -> 0
  lut <- c(2L, NA_integer_, 1L, 0L)
  G <- matrix(lut[codes + 1L], 4 * nbytes, length(snps))
  G[seq_len(n), , drop = FALSE]
}

#' Harmonise summary-statistic alleles to a reference
#'
#' Matches by chromosome and position and returns the sign that aligns each
#' Z-score to the reference A1 allele. Ambiguous (A/T, C/G) SNPs and allele
#' mismatches are dropped.
#'
#' @param chr,pos,effect,other Summary statistics: chromosome, position,
#'   effect allele (e.g. Pan-UKB \code{alt}) and other allele (\code{ref}).
#' @param ref A data frame with \code{chr}, \code{pos}, \code{a1}, \code{a2}
#'   (e.g. from [read_bim()]).
#' @param drop_ambiguous Drop strand-ambiguous SNPs.
#' @return Data frame with \code{i} (row in the summary statistics),
#'   \code{j} (row in \code{ref}) and \code{sign} (+1 / -1).
#' @export
harmonise_alleles <- function(chr, pos, effect, other, ref,
                              drop_ambiguous = TRUE) {
  key_s <- paste(chr, pos, sep = ":")
  key_r <- paste(ref$chr, ref$pos, sep = ":")
  j <- match(key_s, key_r)
  i <- which(!is.na(j))
  j <- j[i]
  e <- toupper(effect[i]); o <- toupper(other[i])
  a1 <- toupper(ref$a1[j]); a2 <- toupper(ref$a2[j])
  sign <- ifelse(e == a1 & o == a2, 1, ifelse(e == a2 & o == a1, -1, NA))
  if (drop_ambiguous) {
    amb <- paste0(e, o) %in% c("AT", "TA", "CG", "GC")
    sign[amb] <- NA
  }
  keep <- !is.na(sign)
  data.frame(i = i[keep], j = j[keep], sign = sign[keep])
}

#' LD reference from per-chromosome PLINK files
#'
#' Builds an \code{lrcpq_ld} object whose \code{get()} computes LD on demand
#' from the reference genotypes (e.g. 1000 Genomes EUR, see
#' \code{inst/scripts/build_1000g_eur_hm3.sh}). Genotypes of one chromosome
#' are cached at a time.
#'
#' @param prefixes Named character vector of PLINK prefixes, names are
#'   chromosomes (e.g. \code{c("22" = "eur_hm3_chr22")}).
#' @param snps Optional data frame of the SNPs to use (columns \code{chr},
#'   \code{pos}), in genomic order; defaults to all reference SNPs. Only
#'   reference SNPs matching \code{snps} are kept.
#' @param blocks LD blocks (default [ld_blocks_eur()]).
#' @param maf_min Drop reference SNPs with lower minor allele frequency.
#' @return An \code{lrcpq_ld} object; \code{$snp} holds \code{chr},
#'   \code{pos}, \code{rsid}, \code{a1}, \code{a2}, \code{maf}, \code{block}.
#' @export
ld_from_plink <- function(prefixes, snps = NULL, blocks = ld_blocks_eur(),
                          maf_min = 0.01) {
  bims <- lapply(names(prefixes), function(cn) {
    b <- read_bim(prefixes[[cn]])
    b$file <- prefixes[[cn]]
    b$row <- seq_len(nrow(b))
    b
  })
  bim <- do.call(rbind, bims)
  if (!is.null(snps)) {
    keep <- paste(bim$chr, bim$pos) %in% paste(snps$chr, snps$pos)
    bim <- bim[keep, ]
  }
  # MAF filter, computed chromosome by chromosome
  bim$maf <- NA_real_
  for (f in unique(bim$file)) {
    sel <- which(bim$file == f)
    G <- read_bed(f, bim$row[sel])
    p <- colMeans(G, na.rm = TRUE) / 2
    bim$maf[sel] <- pmin(p, 1 - p)
  }
  bim <- bim[bim$maf >= maf_min, ]
  bim <- bim[order(bim$chr, bim$pos), ]
  bim$block <- assign_blocks(bim$chr, bim$pos, blocks)
  N <- length(readLines(paste0(bim$file[1], ".fam")))
  cache <- new.env()
  cache$file <- ""
  get <- function(idx) {
    files <- unique(bim$file[idx])
    Gs <- NULL
    for (f in files) {
      if (cache$file != f) {
        sel <- which(bim$file == f)
        G <- read_bed(f, bim$row[sel])
        G <- scale(G)
        G[is.na(G)] <- 0
        cache$G <- G
        cache$map <- sel
        cache$file <- f
      }
      cols <- match(idx[bim$file[idx] == f], cache$map)
      Gs <- cbind(Gs, cache$G[, cols, drop = FALSE])
    }
    R <- crossprod(Gs) / (N - 1)
    diag(R) <- 1
    R
  }
  snp <- bim[, c("chr", "pos", "rsid", "a1", "a2", "maf", "block")]
  rownames(snp) <- NULL
  new_ld(get, nrow(bim), bim$block, snp, N)
}

#' LD reference from a directory of overlapping LD windows
#'
#' Loads the window format produced by the LRCQ analysis scripts from the
#' UK Biobank in-sample LD release (PolyFun, \code{broad-alkesgroup-ukbb-ld}):
#' one \code{<win>.snps.tsv} (with a \code{chr:pos:ref:alt} id column or
#' \code{chr}/\code{pos} columns) and one dense float32 \code{<win>.R.f32}
#' per window, windows named \code{chr<c>_<start>_<end>}. Each window's core
#' is its middle \code{step} bp (edge windows extend their core to the
#' chromosome ends), so cores tile the genome.
#'
#' @param dir Directory holding the window files.
#' @param step Window step in bp (core length).
#' @param id_col Name of the SNP id column in the \code{.snps.tsv} files
#'   (\code{chr:pos:ref:alt}; the LRCQ analysis files use \code{ID}, with
#'   \code{CHR} and \code{BP} columns and R already sign-aligned to ID).
#' @param N_ref Reference panel size used for the r-squared bias correction
#'   (UKB British LD: about 337,000).
#' @return List with \code{ld} (an \code{lrcpq_ld} over the union of SNPs,
#'   whose \code{get()} reads the window covering the request) and
#'   \code{windows} (pass to \code{lrcq(windows = )}).
#' @export
ld_from_windows <- function(dir, step = 2e6, id_col = "ID", N_ref = 337000) {
  files <- list.files(dir, pattern = "\\.snps\\.tsv$", full.names = TRUE)
  if (!length(files)) stop("no *.snps.tsv files in ", dir)
  nm <- sub("\\.snps\\.tsv$", "", basename(files))
  parts <- regmatches(nm, regexec("chr([0-9]+)_([0-9]+)_([0-9]+)", nm))
  meta <- data.frame(name = nm, file = files,
                     chr = as.integer(vapply(parts, `[`, "", 2)),
                     start = as.numeric(vapply(parts, `[`, "", 3)),
                     end = as.numeric(vapply(parts, `[`, "", 4)))
  meta <- meta[order(meta$chr, meta$start), ]
  tabs <- lapply(meta$file, function(f) {
    t <- utils::read.delim(f, stringsAsFactors = FALSE)
    if (!"chr" %in% names(t) && "CHR" %in% names(t)) t$chr <- t$CHR
    if (!"pos" %in% names(t) && "BP" %in% names(t)) t$pos <- t$BP
    if (!id_col %in% names(t)) {
      t[[id_col]] <- paste(t$chr, t$pos, sep = ":")
    }
    if (!"pos" %in% names(t)) {
      t$chr <- as.integer(sub(":.*", "", t[[id_col]]))
      t$pos <- as.integer(sub("^[^:]*:([0-9]+).*", "\\1", t[[id_col]]))
    }
    t
  })
  all_snp <- unique(do.call(rbind, lapply(tabs, function(t) t[, c(id_col, "chr", "pos")])))
  all_snp <- all_snp[order(all_snp$chr, all_snp$pos), ]
  rownames(all_snp) <- NULL
  ids <- all_snp[[id_col]]
  flank <- (meta$end - meta$start - step) / 2
  windows <- vector("list", nrow(meta))
  for (k in seq_len(nrow(meta))) {
    first <- k == 1 || meta$chr[k - 1] != meta$chr[k]
    last <- k == nrow(meta) || meta$chr[k + 1] != meta$chr[k]
    lo <- if (first) -Inf else meta$start[k] + flank[k]
    hi <- if (last) Inf else meta$end[k] - flank[k]
    idx <- match(tabs[[k]][[id_col]], ids)
    pos <- all_snp$pos[idx]
    windows[[k]] <- list(idx = idx, core = which(pos >= lo & pos < hi),
                         file = sub("\\.snps\\.tsv$", ".R.f32", meta$file[k]))
  }
  cache <- new.env(); cache$k <- 0L
  get <- function(idx) {
    k <- which(vapply(windows, function(w) all(idx %in% w$idx), TRUE))[1]
    if (is.na(k)) stop("no single LD window covers the requested SNPs")
    if (cache$k != k) {
      p <- length(windows[[k]]$idx)
      con <- file(windows[[k]]$file, "rb")
      v <- readBin(con, "numeric", n = p * p, size = 4)
      close(con)
      cache$R <- matrix(v, p, p)
      cache$k <- k
    }
    j <- match(idx, windows[[k]]$idx)
    cache$R[j, j, drop = FALSE]
  }
  block <- cumsum(c(TRUE, diff(all_snp$chr) != 0))
  list(ld = new_ld(get, nrow(all_snp), block, all_snp, N = N_ref),
       windows = windows)
}
