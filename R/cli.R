#' Command-line interface
#'
#' Entry point of the \code{lrcpq} command-line tool
#' (\code{Rscript $(Rscript -e 'cat(system.file("cli/lrcpq", package = "lrcpq"))') <command> ...}),
#' modelled on \code{ldsc.py}. Commands:
#' \describe{
#'   \item{\code{munge}}{\code{--sumstats FILE --ref-ld PREFIX@ --out FILE [--N n]}:
#'     one GWAS file to a harmonised table (chr, pos, a1, a2, z, n) aligned to
#'     the LD reference A1 allele ([munge_sumstats()]).}
#'   \item{\code{lrcq}}{\code{--sumstats-list FILE --ref-ld PREFIX@ --chr 1-22 --M m
#'     [--gcov FILE --intercept FILE | --h2 FILE --rg FILE --intercept FILE]
#'     [--method ols|wls] [--skew] [--tag-r2 0.5] --out PREFIX}: stages 1-3 with the
#'     traits listed one munged file per line (name<TAB>path, or a path whose
#'     basename names the trait). Without stage 1-2 files they are estimated by
#'     LD score regression, which needs genome-wide data. Writes
#'     \code{PREFIX.lrcq.tsv.gz}, \code{PREFIX.stage12.rds}, \code{PREFIX.input.rds}.}
#'   \item{\code{lrcp}}{\code{--lrcq PREFIX [--threshold 10] [--z-min 2]
#'     [--min-distance 5e6] [--class-null] [--n-boot 1000] --out PREFIX}: stage 4
#'     on the candidates of an \code{lrcq} run ([run_lrcp()]); writes
#'     \code{PREFIX.lrcp.tsv.gz}.}
#' }
#' \code{--ref-ld} takes a PLINK prefix with \code{@} standing for the
#' chromosome, as in ldsc. Every command writes \code{OUT.log} with the call,
#' package version, input counts and timings.
#' @param args Character vector of arguments (\code{commandArgs(TRUE)}).
#' @return Invisibly, the output written.
#' @export
lrcpq_cli <- function(args = commandArgs(TRUE)) {
  if (!length(args) || args[1] %in% c("-h", "--help", "help")) {
    cat("usage: lrcpq <munge|lrcq|lrcp> [options]; see ?lrcpq_cli\n"); return(invisible(NULL))
  }
  cmd <- args[1]; opt <- cli_parse(args[-1])
  out <- opt$out; if (is.null(out)) stop("--out is required")
  logf <- paste0(if (cmd == "munge") sub("\\.gz$", "", out) else out, ".log")
  log <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S "), ...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
  cat("", file = logf)
  log("lrcpq ", as.character(utils::packageVersion("lrcpq")), " ", cmd, " ", paste(args[-1], collapse = " "))
  t0 <- Sys.time()
  res <- switch(cmd,
    munge = cli_munge(opt, log),
    lrcq = cli_lrcq(opt, log),
    lrcp = cli_lrcp(opt, log),
    stop("unknown command: ", cmd))
  log("done in ", format(round(Sys.time() - t0, 1)))
  invisible(res)
}

cli_parse <- function(a) {
  opt <- list(); i <- 1
  while (i <= length(a)) {
    if (!startsWith(a[i], "--")) stop("unexpected argument: ", a[i])
    key <- gsub("-", "_", substring(a[i], 3))
    if (i < length(a) && !startsWith(a[i + 1], "--")) { opt[[key]] <- a[i + 1]; i <- i + 2 }
    else { opt[[key]] <- TRUE; i <- i + 1 }
  }
  opt
}

cli_chrs <- function(x) {
  if (is.null(x)) return(as.character(1:22))
  unlist(lapply(strsplit(x, ",")[[1]], function(p)
    if (grepl("-", p)) { r <- as.integer(strsplit(p, "-")[[1]]); as.character(r[1]:r[2]) } else p))
}

cli_ref <- function(opt) {
  if (is.null(opt$ref_ld)) stop("--ref-ld PREFIX@ is required")
  chrs <- cli_chrs(opt$chr)
  pre <- stats::setNames(vapply(chrs, function(c) sub("@", c, opt$ref_ld, fixed = TRUE), ""), chrs)
  pre[file.exists(paste0(pre, ".bim"))]
}

cli_munge <- function(opt, log) {
  pre <- cli_ref(opt)
  ref <- do.call(rbind, lapply(pre, read_bim))
  log("reference SNPs: ", nrow(ref), " on ", length(pre), " chromosome(s)")
  m <- munge_sumstats(opt$sumstats, ref, n = if (is.null(opt$N)) NULL else as.numeric(opt$N),
                      info_min = as.numeric(if (is.null(opt$info_min)) 0.9 else opt$info_min),
                      maf_min = as.numeric(if (is.null(opt$maf_min)) 0.01 else opt$maf_min))
  lg <- attr(m, "log"); log(paste(names(lg), lg, sep = ": ", collapse = "; "))
  out <- data.frame(chr = ref$chr[m$j], pos = ref$pos[m$j], rsid = ref$rsid[m$j],
                    a1 = ref$a1[m$j], a2 = ref$a2[m$j], z = signif(m$z, 6), n = m$n)
  con <- if (grepl("\\.gz$", opt$out)) gzfile(opt$out, "w") else file(opt$out, "w")
  utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
  log("wrote ", nrow(out), " SNPs to ", opt$out)
  invisible(out)
}

cli_read_matrix <- function(f, ids) {
  x <- utils::read.delim(f, row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  if (!all(ids %in% rownames(x))) stop(f, " lacks traits: ", paste(setdiff(ids, rownames(x)), collapse = ", "))
  x[ids, ids]
}

cli_lrcq <- function(opt, log) {
  lst <- utils::read.delim(opt$sumstats_list, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(lst) == 1) lst <- data.frame(V1 = sub("\\..*$", "", basename(lst$V1)), V2 = lst$V1)
  ids <- lst$V1
  ss <- lapply(lst$V2, function(f) utils::read.delim(f, stringsAsFactors = FALSE))
  key <- function(d) paste(d$chr, d$pos, sep = ":")
  common <- Reduce(intersect, lapply(ss, key))
  log(length(ids), " traits; ", length(common), " SNPs present in all")
  snps <- do.call(rbind, strsplit(common, ":")); snps <- data.frame(chr = as.integer(snps[, 1]), pos = as.integer(snps[, 2]))
  pre <- cli_ref(opt)
  ld <- ld_from_plink(pre, snps = snps)
  k <- paste(ld$snp$chr, ld$snp$pos, sep = ":")
  Z <- sapply(ss, function(d) { i <- match(k, key(d))
    if (any(d$a1[i] != ld$snp$a1)) stop("munged files must be aligned to the same reference (--ref-ld)")
    d$z[i] })
  colnames(Z) <- ids
  n <- vapply(ss, function(d) stats::median(d$n), 0)
  log("LD reference: ", ld$m, " SNPs in ", length(unique(ld$block_id)), " blocks; N_ref = ", ld$N)
  M <- if (is.null(opt$M)) { log("--M not given: using the number of reference SNPs"); ld$m } else as.numeric(opt$M)
  gcov <- intercept <- NULL
  if (!is.null(opt$intercept)) intercept <- cli_read_matrix(opt$intercept, ids)
  if (!is.null(opt$gcov)) gcov <- cli_read_matrix(opt$gcov, ids)
  if (!is.null(opt$h2) && !is.null(opt$rg)) {
    h2 <- utils::read.delim(opt$h2, header = FALSE, row.names = 1)[ids, 1]
    gcov <- gcov_from_rg(h2, cli_read_matrix(opt$rg, ids))
  }
  ldscore <- NULL
  if (is.null(gcov)) { log("estimating stage 1-2 by LD score regression"); ldscore <- ld_scores(ld) }
  method <- if (is.null(opt$method)) "wls" else opt$method
  skew <- isTRUE(opt$skew)
  if (skew && method != "ols") { log("--skew needs OLS: using --method ols"); method <- "ols" }
  tag_r2 <- if (is.null(opt$tag_r2)) 0.5 else as.numeric(opt$tag_r2)
  fit <- run_lrcq(Z, ld, n, M, gcov = gcov, intercept = intercept, ldscore = ldscore,
                  tag_r2 = tag_r2, method = method, skew = skew)
  sc <- check_scale(Z, ld, n, diag(fit$stage12$gcov), M, diag(as.matrix(fit$stage12$intercept)))
  log("check_scale median ", round(stats::median(sc), 3), " (near 1 when n, h2 and M are on one scale)")
  w <- fit$w; w <- w[!is.na(w$w_raw), ]
  con <- gzfile(paste0(opt$out, ".lrcq.tsv.gz"), "w")
  utils::write.table(w, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
  saveRDS(fit$stage12, paste0(opt$out, ".stage12.rds"))
  saveRDS(list(Z = Z, n = n, M = M, ref_ld = pre, snps = snps, fit = fit), paste0(opt$out, ".input.rds"))
  log("wrote ", nrow(w), " tag estimates to ", opt$out, ".lrcq.tsv.gz")
  invisible(fit)
}

cli_lrcp <- function(opt, log) {
  if (is.null(opt$lrcq)) stop("--lrcq PREFIX (from the lrcq command) is required")
  inp <- readRDS(paste0(opt$lrcq, ".input.rds"))
  ld <- ld_from_plink(inp$ref_ld, snps = inp$snps)
  num <- function(x, d) if (is.null(x)) d else as.numeric(x)
  pairs <- run_lrcp(inp$Z, ld, inp$n, inp$M, inp$fit$stage12$gcov, inp$fit$stage12$intercept,
                    inp$fit, threshold = num(opt$threshold, 10), z_min = num(opt$z_min, -Inf),
                    min_distance = num(opt$min_distance, 5e6), class_null = isTRUE(opt$class_null),
                    n_boot = num(opt$n_boot, 1000))
  if (nrow(pairs)) {
    pairs$rsid1 <- ld$snp$rsid[pairs$snp1]; pairs$rsid2 <- ld$snp$rsid[pairs$snp2]
    pairs <- pairs[order(pairs$p), ]
  }
  con <- gzfile(paste0(opt$out, ".lrcp.tsv.gz"), "w")
  utils::write.table(pairs, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
  log("wrote ", nrow(pairs), " pairs to ", opt$out, ".lrcp.tsv.gz")
  invisible(pairs)
}
