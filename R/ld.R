# LD reference objects, sliding windows and LD-derived design quantities.

#' LD reference object
#'
#' An \code{lrcpq_ld} object stores SNP metadata, a partition of SNPs into
#' approximately independent LD blocks, and a function returning the LD
#' correlation matrix of any set of SNP indices. Constructors:
#' [ld_from_matrix()], [ld_from_blocks()], [ld_from_genotypes()].
#'
#' @param get Function of an integer index vector returning that LD matrix.
#' @param m Number of SNPs.
#' @param block_id Integer block label per SNP (consecutive).
#' @param snp Optional SNP metadata data frame.
#' @param N Reference panel sample size (for r-squared bias correction;
#'   \code{Inf} for population LD).
#' @export
new_ld <- function(get, m, block_id, snp = NULL, N = Inf) {
  stopifnot(length(block_id) == m, !is.unsorted(block_id))
  structure(list(get = get, m = m, block_id = block_id, snp = snp, N = N),
            class = "lrcpq_ld")
}

#' @export
print.lrcpq_ld <- function(x, ...) {
  cat("<lrcpq_ld>", x$m, "SNPs in", length(unique(x$block_id)),
      "blocks; reference N =", x$N, "\n")
  invisible(x)
}

#' LD reference from a full LD matrix
#' @param R m-by-m LD correlation matrix.
#' @param block_sizes Sizes of consecutive LD blocks (default: one block).
#' @inheritParams new_ld
#' @export
ld_from_matrix <- function(R, block_sizes = nrow(R), snp = NULL, N = Inf) {
  stopifnot(sum(block_sizes) == nrow(R))
  new_ld(function(idx) as.matrix(R[idx, idx, drop = FALSE]), nrow(R),
         rep(seq_along(block_sizes), block_sizes), snp, N)
}

#' LD reference from a list of independent LD blocks
#' @param blocks List of LD block matrices (zero LD between blocks).
#' @inheritParams new_ld
#' @export
ld_from_blocks <- function(blocks, snp = NULL, N = Inf) {
  sizes <- vapply(blocks, nrow, 1L)
  R <- Matrix::bdiag(blocks)
  new_ld(function(idx) as.matrix(R[idx, idx, drop = FALSE]), sum(sizes),
         rep(seq_along(sizes), sizes), snp, N)
}

#' LD reference computed on the fly from a genotype matrix
#'
#' @param G Individuals-by-SNPs genotype dosage matrix (0/1/2), SNPs in
#'   genomic order.
#' @param block_id Block label per SNP (see [assign_blocks()]).
#' @inheritParams new_ld
#' @export
ld_from_genotypes <- function(G, block_id, snp = NULL) {
  Gs <- scale(G)
  Gs[is.na(Gs)] <- 0
  N <- nrow(G)
  new_ld(function(idx) {
    R <- crossprod(Gs[, idx, drop = FALSE]) / (N - 1)
    diag(R) <- 1
    R
  }, ncol(G), block_id, snp, N)
}

#' Assign SNPs to LD blocks from block boundaries
#'
#' @param chr,pos SNP chromosome and base-pair position (sorted).
#' @param blocks Data frame with columns \code{chr}, \code{start}, \code{stop}
#'   (e.g. Berisa and Pickrell 2016 EUR blocks).
#' @return Integer block id per SNP, consecutive.
#' @export
assign_blocks <- function(chr, pos, blocks) {
  id <- rep(NA_integer_, length(pos))
  for (b in seq_len(nrow(blocks))) {
    sel <- chr == blocks$chr[b] & pos >= blocks$start[b] & pos < blocks$stop[b]
    id[sel] <- b
  }
  if (anyNA(id)) {
    # SNPs outside listed blocks join the nearest preceding block
    for (i in which(is.na(id))) id[i] <- if (i > 1) id[i - 1] else 1L
  }
  as.integer(factor(id, levels = unique(id)))
}

#' Sliding windows of consecutive LD blocks
#'
#' Chapter 6 strategy: each window spans \code{size} consecutive blocks
#' (two at chromosome ends), slides one block at a time, and only the
#' estimates for the middle (core) block are kept.
#'
#' @param block_id Block label per SNP.
#' @param size Number of blocks per window (odd).
#' @param chr Optional chromosome per SNP; windows never cross chromosomes.
#' @return List of windows, each with \code{idx} (all SNPs) and \code{core}
#'   (positions within \code{idx} of the core block SNPs).
#' @export
make_windows <- function(block_id, size = 3, chr = NULL) {
  stopifnot(size %% 2 == 1)
  half <- (size - 1) / 2
  blocks <- unique(block_id)
  bchr <- if (is.null(chr)) rep(1, length(blocks)) else chr[match(blocks, block_id)]
  lapply(seq_along(blocks), function(b) {
    nb <- seq(max(1, b - half), min(length(blocks), b + half))
    nb <- nb[bchr[nb] == bchr[b]]
    idx <- which(block_id %in% blocks[nb])
    list(idx = idx, core = which(block_id[idx] == blocks[b]), block = blocks[b])
  })
}

#' Squared LD with small-sample bias correction
#'
#' \eqn{\tilde r^2 = r^2 - (1 - r^2)/(N - 2)} off the diagonal (Bulik-Sullivan
#' et al. 2015), diagonal kept at 1.
#' @param R LD correlation matrix.
#' @param N Reference sample size (\code{Inf}: no correction).
#' @export
ld_r2 <- function(R, N = Inf) {
  D <- as.matrix(R)^2
  if (is.finite(N)) {
    dg <- diag(D)
    D <- D - (1 - D) / (N - 2)
    diag(D) <- dg
  }
  D
}

#' LD scores within an LD reference
#'
#' \eqn{l_k = \sum_s \tilde r^2_{ks}}, summed within each sliding window
#' (blocks plus neighbours) so edge SNPs see their cross-block LD.
#' @param ld An \code{lrcpq_ld} object.
#' @param size Window size in blocks.
#' @export
ld_scores <- function(ld, size = 3) {
  out <- numeric(ld$m)
  for (win in make_windows(ld$block_id, size)) {
    D <- ld_r2(ld$get(win$idx), ld$N)
    out[win$idx[win$core]] <- rowSums(D[win$core, , drop = FALSE])
  }
  out
}

#' Prune SNPs in near-perfect LD
#'
#' Greedy pruning keeping the first SNP of any pair with \eqn{r^2 >}
#' \code{r2_max}. Near-duplicate SNPs make \eqn{D = R \circ R} singular;
#' their enrichment is not separately identifiable.
#' @param R LD matrix.
#' @param r2_max Threshold.
#' @return Indices of kept SNPs.
#' @export
prune_ld <- function(R, r2_max = 0.9) {
  R2 <- as.matrix(R)^2
  keep <- logical(nrow(R2))
  dropped <- logical(nrow(R2))
  for (i in seq_len(nrow(R2))) {
    if (dropped[i]) next
    keep[i] <- TRUE
    dropped <- dropped | (R2[i, ] > r2_max)
  }
  which(keep)
}
