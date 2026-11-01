# lrcpq

Phenome-wide genetic-effect enrichment (LRCQ) and genetic-effect correlation (LRCP) from GWAS summary statistics and an LD reference. The package implements this unified matrix-normal framework:

- **LRCQ**: per-SNP phenome-wide genetic-effect enrichment `w` from squared Z-scores of many GWAS.
- **LRCP**: genetic-effect correlation `Rβ` between SNPs (or genes) from Z-score cross-products.

Model: `Y = XB + E`, `B ~ MN(0, U, V)` with `U = W^{1/2} Rβ W^{1/2}` and `V = H_g / M`; GWAS Z-scores follow `Z = R B diag(√n) + E`, `E ~ MN(0, R, C)`, where `C` holds the inter-GWAS intercepts (`r_ab · o_ab`).

## Install

```r
# install.packages("remotes")
remotes::install_github("unfated/LRCPQ", build_vignettes = TRUE)
```

The repository is private until the papers are published. Until then, install from a clone with `R CMD INSTALL .`

## Quick start

`vignette("getting-started", package = "lrcpq")` runs the whole pipeline (munge, then LRCQ, then LRCP) on the bundled chr22 example in a few minutes. From the shell, use the `ldsc.py`-style command-line tool:

```sh
LRCPQ=$(Rscript -e 'cat(system.file("cli","lrcpq",package="lrcpq"))')
$LRCPQ munge --sumstats trait.sumstats.gz --ref-ld ref/chr@ --chr 1-22 --out trait.munged.gz
$LRCPQ lrcq  --sumstats-list traits.txt --ref-ld ref/chr@ --chr 1-22 --M 1094844 \
             --h2 h2.tsv --rg rg.tsv --intercept intercept.tsv --method ols --skew --out run
$LRCPQ lrcp  --lrcq run --out run
```

Each command writes `OUT.log`. See `?lrcpq_cli` for every flag.

## Input formats

| Input | Format |
|---|---|
| GWAS summary statistics | Any whitespace- or tab-delimited file, optionally gzipped, with SNP or CHR+BP, both alleles, and one of Z, BETA+SE, OR+SE or signed P, plus N or N_cases+N_controls. Column names follow LDSC's synonyms; override them with `cols`. |
| LD reference | PLINK `.bed/.bim/.fam` per chromosome (`@` = chromosome). The prebuilt 1000G EUR HapMap3 reference (503 samples, GRCh37) is attached to the [v0.2.0 release](https://github.com/unfated/LRCPQ/releases/tag/v0.2.0) and can be rebuilt with `inst/scripts/build_1000g_eur_hm3.sh`. Alternatively, use UK Biobank in-sample LD windows from the public PolyFun/Broad bucket (`ld_from_windows()`). |
| `--sumstats-list` | Two columns: trait name and munged file path. |
| `--h2` | Two columns, no header: trait name and SNP heritability (the same LDSC fit as N). |
| `--rg`, `--intercept`, `--gcov` | Square tab-delimited matrices with trait names as row and column names. |

## Outputs

- `*.lrcq.tsv.gz` has one row per tag SNP, with these columns:
  - `w_raw`, the clump-level enrichment, and its `se`;
  - `w_null1`, the expected estimate with no enrichment;
  - with `--skew`: `k3`, `p_gt0` (H0: w = 0) and `p_gt1` (H0: no enrichment).
- `*.lrcp.tsv.gz` has one row per distal tag pair: `rho`, `C`, `se`, `z`, bootstrap `p`, `z_gls`, `q_BH`, `q_BY`, and the rsids.

Speed: use an optimised BLAS. With OpenBLAS (`apt-get install libopenblas0-pthread`) instead of R's reference BLAS, per-tag LRCQ on chromosome 22 (15,850 HM3 SNPs, 100 traits, 1000G LD) drops from 12 minutes to 42 seconds, roughly 50 minutes genome-wide.

## Four-stage estimator

| Stage | Parameter | Function |
|---|---|---|
| 1 | `h²`, `h_ab` | `ldsc_h2()`, `ldsc_gcov()`, `ldsc_matrix()` |
| 2 | intercepts `C` | `ldsc_matrix()` (cross-trait intercepts), `make_intercept()` for known overlap |
| 3 | `w` (LRCQ) | `lrcq()` / `lrcq_window()`; `method = "ols" / "wls" / "irls" / "equal"`; `rectify_w()` methods A, B, C |
| 4 | `Rβ` (LRCP) | `lrcp_distal()` (GLS/OLS/WLS, `se_type = "null"/"plugin"`), `lrcp_mle()`, `lrcp_local()`; gene level `lrcp_gene()` |

Stage 3 uses the fast OLS/WLS solution: stacking q traits gives the design `s ⊗ D` (`D = R∘R`), so the regression collapses to one m-by-m solve per LD window, `(D Ω D) w = D u`.

## Real data

| Task | Functions |
|---|---|
| Pan-UKB phenotypes and Z-scores (tabix + HTTP range requests, no full download) | `panukb_manifest()`, `panukb_select()`, `panukb_z()` |
| LD reference | `ld_from_plink()` (1000G EUR, built by `inst/scripts/build_1000g_eur_hm3.sh`), `ld_from_windows()` (UKB in-sample LD windows), `ld_blocks_eur()`, `harmonise_alleles()` |
| Munging | `munge_sumstats()`, `build_z_matrix()`, command line `lrcpq munge` |
| Pipelines | `run_lrcq()` (stages 1-3), `run_lrcp()` (screen on `w`, distal block pairs, BH/BY q-values), `gcov_from_rg()`, `prune_traits()`, `check_scale()` (n, h2 and M on one scale), `lrcp_screen()` (genome-wide screen with simulated-null calibration), `p_skew()` (skew-corrected LRCQ p) |
| Aggregation | `lrcq_category()`, `lrcq_annot_regression()`, `trait_weights()` |
| Power | `power_lrcq()`, `power_lrcp()`, `detectable_effect()` |

Vignettes: `vignette("getting-started")`, `vignette("simulation")` and `vignette("real-data")` (package `lrcpq`).

## Simulation in five lines

```r
library(lrcpq)
blocks <- make_ld_ar1(rep(50, 10), rho = 0.5)          # toy LD, or real LD blocks
ld     <- ld_from_blocks(blocks)
sim    <- simulate_lrcpq(blocks, n = rep(2e4, 30), prop_nonzero = 0.1)
fit    <- lrcq(sim$Z, ld, sim$n, sim$gcov, M = ld$m, intercept = sim$intercept)
cor(fit$w_raw, sim$w)
```

## FAQ

**What is M?** M is the number of SNPs that heritability is spread over in the LDSC fit behind h2: 1,094,844 for the HapMap3 LD-score set. It is not LDSC's `M_5_50`. `check_scale()` checks that n, h2 and M agree. Its ratio should be about 1 genome-wide.

**Which N?** Use total N from the same fit as h2, not effective N.

**Where do h2, rg and intercepts come from?** From genome-wide LD score regression: `ldsc_matrix()`, published LDSC output, or `panukb_select()` for Pan-UKB. Use the genome-wide rg (`gcov_from_rg()`). A diagonal gcov makes SEs far too small for correlated traits, so the package warns about it. The intercept matrix must come from one joint fit, so that it is positive semi-definite.

**Why skew-corrected p-values?** Phenome-wide chi-square sums are right-skewed. The normal reference is anti-conservative in the upper tail, so `p_skew()` uses a shifted gamma with the estimate's third cumulant.

**Why are LRCP bootstrap p-values conservative for strong pairs?** The statistic is self-normalised and bounded. Its tails are lighter than normal, so the normal p overstates significance less often than it understates it. Report `rho` descriptively and use `lrcp_screen()` with the simulated null for genome-wide FDR.

**Why no Tukey range or sample-variance test across traits?** Both assume independent traits. Phenome-wide traits share genetic correlation and samples, so these tests are invalid here. The trait-weighted chi-square of LRCQ accounts for both.

**How fast is it?** Use OpenBLAS (see above). Genome-wide LRCQ for about 450 traits runs in hours, not days. Split long runs by chromosome.

## Unsettled maths, kept switchable

- `pair_sum = "once"` (default, confirmed by the theory supplement) counts each SNP pair once in the `Rβ` cross terms; `"chapter"` reproduces the thesis equations, which count each pair twice.
- `lrcq(..., cross_terms = FALSE)` is the chi-square-only LRCQ of the abstract; `TRUE` adds the cross terms for a given `Rβ`.

Intercepts are taken from stage 2 and held fixed: with a free per-trait intercept inside a window, `w` is not identifiable.

- `se_rectify` in `lrcq_window()` sets the plug-in for the LRCQ model SE. The default is `"bgA"`, background-filled rule A as in the LRCQ paper (theory supplement S3.10a). `"A"` was the 0.1 default; it drops the polygenic background and makes SEs 10-25% too small under dense enrichment.

## Citation

`citation("lrcpq")`. Licence: GPL (>= 3).
