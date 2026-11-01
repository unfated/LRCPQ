# Pan-UK Biobank summary statistics: manifest, tabix index, range-request
# region fetching (adapted from the LRCP analysis thread's approach: parse the
# .tbi index and issue one HTTP Range request per region, instead of remote
# tabix through htslib, which over-reads and truncates on S3).

#' Read the Pan-UKB phenotype manifest
#'
#' @param path Local path or URL of \code{phenotype_manifest.tsv.bgz}
#'   (default: the public S3 copy).
#' @return Data frame.
#' @export
panukb_manifest <- function(path = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz") {
  if (grepl("^https?://", path)) {
    tmp <- tempfile(fileext = ".tsv.gz")
    http_get(path, tmp)
    path <- tmp
  }
  utils::read.delim(gzfile(path), stringsAsFactors = FALSE, quote = "")
}

#' Select Pan-UKB phenotypes for a phenome-wide analysis
#'
#' @param manifest From [panukb_manifest()].
#' @param pop Ancestry group whose QC must pass.
#' @param sex Sex stratum (\code{"both_sexes"}).
#' @param independent_only Keep only the maximal independent set.
#' @param min_n Minimum EUR sample size (cases + controls).
#' @return The selected manifest rows with column \code{n}, the total sample
#'   size (cases + controls). Use it as \code{n} with Pan-UKB's observed-scale
#'   h2: s_a = n h2 / M needs n from the same fit as h2, which for Pan-UKB LDSC is
#'   the total N, not an effective N (4 p (1 - p) N for binary traits).
#'   \code{n_eff} is kept as an alias of \code{n} for old scripts; despite its
#'   name it is the total N.
#' @export
panukb_select <- function(manifest, pop = "EUR", sex = "both_sexes",
                          independent_only = FALSE, min_n = 0) {
  qc <- manifest[[paste0("phenotype_qc_", pop)]]
  keep <- !is.na(qc) & qc == "PASS" & manifest$pheno_sex == sex
  if (independent_only) keep <- keep & manifest$in_max_independent_set %in% c(TRUE, "true", "True")
  nc <- suppressWarnings(as.numeric(manifest[[paste0("n_cases_", pop)]]))
  nn <- suppressWarnings(as.numeric(manifest[[paste0("n_controls_", pop)]]))
  n <- nc + ifelse(is.na(nn), 0, nn)
  keep <- keep & !is.na(n) & n >= min_n
  out <- manifest[keep, ]
  out$n <- n[keep]
  out$n_eff <- out$n   # deprecated alias: total N, not an effective N
  out
}

#' HTTP GET (optionally a byte range) with retries via the curl binary
#' @keywords internal
http_get <- function(url, dest, range = NULL, tries = 4) {
  args <- c("-sfL", "--max-time", "600", "-o", shQuote(dest))
  if (!is.null(range)) args <- c(args, "-r", sprintf("%.0f-%.0f", range[1], range[2]))
  for (i in seq_len(tries)) {
    st <- system2("curl", c(args, shQuote(url)))
    if (st == 0) return(invisible(dest))
    Sys.sleep(2^i)
  }
  stop("download failed: ", url)
}

#' Decompress a buffer of whole BGZF blocks
#' @keywords internal
bgzf_inflate <- function(raw) {
  out <- list()
  pos <- 1L
  n <- length(raw)
  while (pos + 17L <= n) {
    if (raw[pos] != as.raw(0x1f) || raw[pos + 1L] != as.raw(0x8b)) stop("bad BGZF block")
    bsize <- readBin(raw[(pos + 16L):(pos + 17L)], "integer", size = 2,
                     signed = FALSE, endian = "little") + 1L
    if (pos + bsize - 1L > n) break        # incomplete trailing block
    blk <- raw[pos:(pos + bsize - 1L)]
    out[[length(out) + 1L]] <- memDecompress(blk, type = "gzip")
    pos <- pos + bsize
  }
  list(data = unlist(out, use.names = FALSE), consumed = pos - 1L)
}

#' Read a tabix (.tbi) index
#'
#' @param path Local path or URL of the \code{.tbi}.
#' @return List with sequence names and, per sequence, the bin chunks and
#'   linear index (virtual file offsets as doubles).
#' @export
read_tbi <- function(path) {
  if (grepl("^https?://", path)) {
    tmp <- tempfile(fileext = ".tbi")
    http_get(path, tmp)
    path <- tmp
  }
  con <- gzfile(path, "rb")
  on.exit(close(con))
  rd_i <- function(k = 1) readBin(con, "integer", k, size = 4, endian = "little")
  rd_u64 <- function(k = 1) {
    v <- readBin(con, "integer", 2 * k, size = 4, endian = "little")
    lo <- v[c(TRUE, FALSE)]; hi <- v[c(FALSE, TRUE)]
    (ifelse(lo < 0, lo + 2^32, lo)) + ifelse(hi < 0, hi + 2^32, hi) * 2^32
  }
  magic <- readBin(con, "raw", 4)
  if (!identical(magic, charToRaw("TBI\001"))) stop("not a tabix index")
  hdr <- rd_i(8)
  n_ref <- hdr[1]
  l_nm <- hdr[8]
  nm <- readBin(con, "raw", l_nm)
  # names are NUL-separated
  nmchr <- rawToChar(replace(nm, nm == as.raw(0), as.raw(0x0a)))
  names <- strsplit(nmchr, "\n")[[1]]
  refs <- vector("list", n_ref)
  for (r in seq_len(n_ref)) {
    n_bin <- rd_i()
    bins <- list()
    for (b in seq_len(n_bin)) {
      bin <- rd_i()
      if (bin < 0) bin <- bin + 2^32
      n_chunk <- rd_i()
      ch <- matrix(rd_u64(2 * n_chunk), ncol = 2, byrow = TRUE)
      bins[[as.character(bin)]] <- ch
    }
    n_intv <- rd_i()
    ioff <- if (n_intv > 0) rd_u64(n_intv) else numeric(0)
    refs[[r]] <- list(bins = bins, ioff = ioff)
  }
  list(names = names, refs = refs, meta_char = hdr[6], skip = hdr[7],
       col_seq = hdr[3], col_beg = hdr[4])
}

#' UCSC binning scheme bins overlapping [beg, end) (0-based)
#' @keywords internal
reg2bins <- function(beg, end) {
  end <- end - 1
  bins <- 0
  for (sh in list(c(1, 26), c(9, 23), c(73, 20), c(585, 17), c(4681, 14))) {
    bins <- c(bins, seq(sh[1] + (beg %/% 2^sh[2]), sh[1] + (end %/% 2^sh[2])))
  }
  bins
}

#' Fetch the rows of a bgzipped, tabix-indexed file overlapping a region
#'
#' Uses the \code{.tbi} bins and linear index to find the compressed byte
#' span, downloads it with one HTTP Range request, inflates the BGZF blocks,
#' and keeps rows inside the region.
#'
#' @param url URL of the \code{.tsv.bgz}.
#' @param tbi Index from [read_tbi()].
#' @param chr Chromosome name as in the file (e.g. \code{"17"}).
#' @param start,end 1-based inclusive region.
#' @param header Optional column names.
#' @return Data frame (character columns unless converted by the caller).
#' @export
fetch_region <- function(url, tbi, chr, start, end, header = NULL) {
  r <- match(as.character(chr), tbi$names)
  if (is.na(r)) return(NULL)
  ref <- tbi$refs[[r]]
  beg0 <- start - 1
  bins <- reg2bins(beg0, end)
  ch <- do.call(rbind, ref$bins[as.character(bins[as.character(bins) %in% names(ref$bins)])])
  if (is.null(ch) || !nrow(ch)) return(NULL)
  li <- min(length(ref$ioff), beg0 %/% 16384 + 1)
  min_off <- if (length(ref$ioff)) ref$ioff[li] else 0
  ch <- ch[ch[, 2] > min_off, , drop = FALSE]
  if (!nrow(ch)) return(NULL)
  vbeg <- max(min(ch[, 1]), min_off)
  vend <- max(ch[, 2])
  cbeg <- floor(vbeg / 65536); ubeg <- vbeg %% 65536
  cend <- floor(vend / 65536)
  tmp <- tempfile()
  on.exit(unlink(tmp))
  http_get(url, tmp, c(cbeg, cend + 65536 + 65535))
  raw <- readBin(tmp, "raw", file.info(tmp)$size)
  txt <- bgzf_inflate(raw)$data
  if (ubeg > 0) txt <- txt[-seq_len(ubeg)]
  lines <- strsplit(rawToChar(txt), "\n", fixed = TRUE)[[1]]
  lines <- lines[!startsWith(lines, "chr\t")]
  if (!length(lines)) return(NULL)
  f <- strsplit(lines, "\t", fixed = TRUE)
  nf <- lengths(f)
  f <- f[nf == stats::median(nf)]
  df <- as.data.frame(do.call(rbind, f), stringsAsFactors = FALSE)
  if (!is.null(header) && length(header) == ncol(df)) names(df) <- header
  sc <- tbi$col_seq; bc <- tbi$col_beg
  pos <- suppressWarnings(as.numeric(df[[bc]]))
  df[df[[sc]] == as.character(chr) & !is.na(pos) & pos >= start & pos <= end, , drop = FALSE]
}

#' Header line of a bgzipped file
#' @param url URL of the \code{.tsv.bgz}.
#' @export
fetch_header <- function(url) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  http_get(url, tmp, c(0, 131071))
  raw <- readBin(tmp, "raw", file.info(tmp)$size)
  txt <- rawToChar(bgzf_inflate(raw)$data)
  strsplit(strsplit(txt, "\n", fixed = TRUE)[[1]][1], "\t", fixed = TRUE)[[1]]
}

#' Pan-UKB Z-scores for a set of SNPs in given regions
#'
#' For each selected phenotype, fetches the regions by HTTP range requests
#' and returns the EUR Z-scores (\code{beta_EUR / se_EUR}) for the requested
#' SNPs, aligned to their \code{chr:pos:ref:alt} ids (Pan-UKB betas are for
#' the \code{alt} allele, GRCh37).
#'
#' @param phenos Rows of the manifest ([panukb_select()]).
#' @param snps Character vector of \code{chr:pos:ref:alt} ids.
#' @param regions Data frame with \code{chr}, \code{start}, \code{end}
#'   (default: one region per chromosome spanning the requested SNPs, padded
#'   by 1 bp).
#' @param base_url Bucket prefix prepended to relative \code{aws_path}s.
#' @param pop Ancestry suffix of the beta/se columns.
#' @param verbose Print progress.
#' @return Matrix of Z-scores (SNPs by phenotypes), \code{NA} where missing.
#' @export
panukb_z <- function(phenos, snps, regions = NULL,
                     base_url = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/",
                     pop = "EUR", verbose = TRUE) {
  sp <- do.call(rbind, strsplit(snps, ":", fixed = TRUE))
  chr <- sp[, 1]; pos <- as.numeric(sp[, 2])
  if (is.null(regions)) {
    regions <- do.call(rbind, lapply(split(pos, chr), function(p) c(min(p), max(p))))
    regions <- data.frame(chr = rownames(regions), start = regions[, 1], end = regions[, 2])
  }
  Z <- matrix(NA_real_, length(snps), nrow(phenos),
              dimnames = list(snps, paste(phenos$trait_type, phenos$phenocode,
                                          phenos$coding, phenos$modifier, sep = "-")))
  for (j in seq_len(nrow(phenos))) {
    url <- phenos$aws_path[j]
    if (!grepl("^https?://", url)) {
      url <- sub("^s3://pan-ukb-us-east-1/sumstats_flat_files/", base_url, url)
      if (!grepl("^https?://", url)) url <- paste0(base_url, phenos$filename[j])
    }
    tbi <- read_tbi(paste0(url, ".tbi"))
    hdr <- fetch_header(url)
    for (k in seq_len(nrow(regions))) {
      d <- fetch_region(url, tbi, regions$chr[k], regions$start[k], regions$end[k], hdr)
      if (is.null(d) || !nrow(d)) next
      id <- paste(d$chr, d$pos, d$ref, d$alt, sep = ":")
      b <- suppressWarnings(as.numeric(d[[paste0("beta_", pop)]]))
      s <- suppressWarnings(as.numeric(d[[paste0("se_", pop)]]))
      hit <- match(snps, id)
      ok <- !is.na(hit)
      Z[ok, j] <- b[hit[ok]] / s[hit[ok]]
    }
    if (verbose) message(sprintf("[%d/%d] %s: %d SNPs", j, nrow(phenos),
                                 colnames(Z)[j], sum(!is.na(Z[, j]))))
  }
  Z
}
