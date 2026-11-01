# lrcpq 0.2.0

* The package moved to its own repository, unfated/LRCPQ.
* New: `munge_sumstats()` and `build_z_matrix()` take raw GWAS files (LDSC column synonyms; Z from BETA/SE, OR/SE or signed P), then filter and align them to the LD reference.
* New: an `ldsc.py`-style command-line tool (`lrcpq_cli()`, `inst/cli/lrcpq`) with `munge`, `lrcq` and `lrcp` commands.
* New: a bundled chr22 example (50 simulated GWAS, 1000G EUR LD), used by the getting-started vignette.
* New: skew-corrected LRCQ p-values (`p_skew()`, `lrcq(..., skew = TRUE)`).
* New: a genome-wide LRCP screen with simulated-null calibration (`lrcp_screen()`).
* New: `check_scale()` checks that n, h2 and M are on one scale.
* Fix: `p_gt1` now tests against `w_null1`, the expected tag estimate with no enrichment. It previously used the clump size.
* Change: LRCQ model SEs (and third cumulants) now use the background-filled rule-A plug-in by default (`se_rectify = "bgA"`), matching the LRCQ paper (theory supplement S3.10a). The old default `"A"` drops the polygenic background, which makes SEs 10-25% too small under dense enrichment. It is still available, along with `"eb"` and `"one"`. Estimates do not change.
* New: the prebuilt 1000G EUR HapMap3 PLINK reference is attached to the GitHub release.
* `panukb_select()` returns total `n` (`n_eff` is kept as an alias).
* `lrcp_gene()` reports bootstrap SDs, and `check_gcov()` warns when the intercept matrix is not positive semi-definite.

# lrcpq 0.1.0

* First version: simulation engine, LD utilities, LD score regression, and the LRCQ and LRCP estimators.
