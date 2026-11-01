# Column-name synonyms (lower case), as in ldsc's munge_sumstats.py
.munge_cols <- list(
  snp = c("snp", "rsid", "rs", "snpid", "marker", "markername", "variant_id", "id"),
  chr = c("chr", "chrom", "chromosome", "#chrom", "#chr"),
  pos = c("pos", "bp", "position", "base_pair_location", "bp_hg19", "pos_hg19"),
  a1 = c("a1", "effect_allele", "alt", "allele1", "ea", "eff_allele", "tested_allele"),
  a2 = c("a2", "other_allele", "ref", "allele2", "nea", "non_effect_allele"),
  z = c("z", "zscore", "z_score", "zstat"),
  beta = c("beta", "b", "effect", "beta_meta", "beta_eur"),
  or = c("or", "odds_ratio"),
  se = c("se", "stderr", "standard_error", "se_meta", "se_eur"),
  p = c("p", "pval", "p_value", "pvalue", "p.value", "pval_meta", "pval_eur"),
  n = c("n", "n_total", "sample_size", "neff_total"),
  n_cas = c("n_cas", "n_cases", "ncase", "n_case"),
  n_con = c("n_con", "n_controls", "ncontrol", "n_control"),
  info = c("info", "imputation_info", "r2", "info_score"),
  maf = c("maf", "frq", "freq", "af", "eaf", "a1_freq", "af_eur", "effect_allele_frequency")
)

#' Read and harmonise one GWAS summary-statistics file
#'
#' The package's counterpart of ldsc's \code{munge_sumstats.py}. Column names
#' are recognised from common synonyms (case-insensitive; override with
#' \code{cols}). The Z-score is taken from a \code{Z} column, else
#' \code{BETA/SE}, \code{log(OR)/SE}, or the p-value with the sign of
#' \code{BETA}/\code{OR}. SNPs are matched to the LD reference by chromosome and
#' position (or by \code{SNP} id when positions are missing), aligned so that Z
#' refers to the reference A1 allele, and filtered on INFO, MAF and
#' strand-ambiguity. Use the reference's genome build (GRCh37 for the bundled
#' 1000G and UKB references).
#'
#' @param x A file path (plain or gzipped, whitespace- or tab-delimited) or a
#'   data frame.
#' @param ref Reference SNP table with \code{chr}, \code{pos}, \code{a1},
#'   \code{a2} and optionally \code{rsid} (e.g. \code{ld$snp}).
#' @param n Sample size if the file has no N column. For binary traits give the
#'   total N (cases + controls), the n that matches observed-scale h2.
#' @param cols Named character vector overriding detected columns, e.g.
#'   \code{c(a1 = "ALT", a2 = "REF", z = "tstat")}.
#' @param info_min,maf_min Filters (applied when the column exists).
#' @param drop_ambiguous Drop A/T and C/G SNPs.
#' @return Data frame with one row per reference SNP found: \code{j} (row in
#'   \code{ref}), \code{z}, \code{n}; attribute \code{"log"} with counts at each
#'   step.
#' @export
munge_sumstats <- function(x, ref, n = NULL, cols = NULL, info_min = 0.9,
                           maf_min = 0.01, drop_ambiguous = TRUE) {
  d <- if (is.character(x)) {
    if (requireNamespace("data.table", quietly = TRUE))
      as.data.frame(data.table::fread(x, data.table = FALSE))
    else utils::read.table(x, header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  } else as.data.frame(x)
  lg <- c(read = nrow(d))
  nm <- tolower(names(d))
  pick <- function(key) {
    if (!is.null(cols) && key %in% names(cols)) return(cols[[key]])
    k <- match(.munge_cols[[key]], nm); k <- k[!is.na(k)]
    if (length(k)) names(d)[k[1]] else NULL
  }
  col <- lapply(stats::setNames(names(.munge_cols), names(.munge_cols)), pick)
  get <- function(key) if (is.null(col[[key]])) NULL else d[[col[[key]]]]
  if (is.null(col$a1) || is.null(col$a2)) stop("effect and other allele columns not found; use cols")
  # Z
  z <- get("z")
  if (is.null(z)) {
    b <- get("beta"); if (is.null(b) && !is.null(get("or"))) b <- log(as.numeric(get("or")))
    se <- get("se")
    if (!is.null(b) && !is.null(se)) z <- as.numeric(b) / as.numeric(se)
    else if (!is.null(b) && !is.null(get("p")))
      z <- sign(as.numeric(b)) * stats::qnorm(as.numeric(get("p")) / 2, lower.tail = FALSE)
    else stop("need Z, BETA/OR and SE, or BETA/OR and P")
  }
  z <- as.numeric(z)
  # N
  nn <- get("n")
  if (is.null(nn) && !is.null(get("n_cas")) && !is.null(get("n_con")))
    nn <- as.numeric(get("n_cas")) + as.numeric(get("n_con"))
  if (is.null(nn)) { if (is.null(n)) stop("no N column: give n"); nn <- rep(n, nrow(d)) }
  nn <- as.numeric(nn)
  keep <- is.finite(z) & is.finite(nn) & nn > 0
  if (!is.null(get("info"))) keep <- keep & !(as.numeric(get("info")) < info_min) %in% TRUE
  if (!is.null(get("maf"))) {
    f <- as.numeric(get("maf")); keep <- keep & !(pmin(f, 1 - f) < maf_min) %in% TRUE
  }
  lg["after_filters"] <- sum(keep)
  # match to reference
  if (!is.null(col$chr) && !is.null(col$pos)) {
    chr <- sub("^chr", "", as.character(get("chr")), ignore.case = TRUE)
    h <- harmonise_alleles(chr, get("pos"), get("a1"), get("a2"), ref, drop_ambiguous)
  } else if (!is.null(col$snp) && !is.null(ref$rsid)) {
    j <- match(get("snp"), ref$rsid)
    h <- harmonise_alleles(ifelse(is.na(j), NA, ref$chr[j]), ifelse(is.na(j), NA, ref$pos[j]),
                           get("a1"), get("a2"), ref, drop_ambiguous)
  } else stop("need CHR and POS columns, or SNP ids with ref$rsid")
  h <- h[keep[h$i], ]
  h <- h[!duplicated(h$j), ]
  lg["matched_aligned"] <- nrow(h)
  out <- data.frame(j = h$j, z = h$sign * z[h$i], n = nn[h$i])
  out <- out[order(out$j), ]
  attr(out, "log") <- lg
  out
}

#' Build the SNP-by-trait Z matrix from munged summary statistics
#'
#' @param munged Named list of [munge_sumstats()] outputs, one per trait.
#' @param m Number of reference SNPs (\code{ld$m}).
#' @param complete Keep only SNPs present in every trait (recommended: the
#'   estimators need a complete matrix).
#' @return List with \code{Z} (m-by-q, \code{NA} where missing),
#'   \code{n} (per-trait median N), \code{keep} (reference rows with complete
#'   data when \code{complete = TRUE}).
#' @export
build_z_matrix <- function(munged, m, complete = TRUE) {
  q <- length(munged)
  Z <- matrix(NA_real_, m, q, dimnames = list(NULL, names(munged)))
  n <- numeric(q)
  for (a in seq_len(q)) {
    Z[munged[[a]]$j, a] <- munged[[a]]$z
    n[a] <- stats::median(munged[[a]]$n)
  }
  keep <- if (complete) which(rowSums(is.na(Z)) == 0) else seq_len(m)
  list(Z = Z, n = stats::setNames(n, names(munged)), keep = keep)
}
