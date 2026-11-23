# RareSV

RareSV tests whether rare structural variants (SVs) that overlap coding sequence are associated with a molecular trait. It can adjust for gene-specific loss-of-function (LoF) carrier status and report single-SV follow-up results.

## Terra scan

The [WDL 1.0 workflow](wdl/sv_skato.wdl) reads chromosome-specific PLINK 2 hard calls with `pgenlibr`. It runs LoF-adjusted SKAT-O for `DEL`, `DUP`, and `DEL+DUP`. The gather task combines the mask p-values per gene with the equal-weight Cauchy method and applies 5% FDR across chromosomes.

See the [Terra setup guide](docs/terra_sv_skato.md) for the manifest format, Docker image, input files, and output tables.

## VCF scan and single-SV follow-up in R

Use the [R helper](R/sv_set_tests.R) with a bgzip-compressed, tabix-indexed VCF. The molecular-trait TSV has one row per gene and one column per participant. The LoF TSV needs `sample_id`, `gene_id`, and `has_lof_variant`.

```r
source("R/sv_set_tests.R")

expression <- read_phenotype_table("expression_zscores.tsv", id_col = "gene_id")
lof <- readr::read_tsv(
  "lof_carriers.tsv",
  col_types = readr::cols(.default = readr::col_character())
)

results <- run_sv_cds_tests(
  phenotype_df = expression,
  gene_ids = unique(expression$gene_id),
  vcf_file = "structural_variants.vcf.gz",
  gtf = "coding_annotation.gtf.gz",
  test = "SKATO",
  sv_types = c("DEL", "DUP"),
  lof_table = lof,
  lof_missing_as_noncarrier = TRUE,
  ignore_gene_version = TRUE,
  follow_up = TRUE,
  keep_carrier_scores = TRUE
)

results$summary
results$single_variant
```

Use `lof_missing_as_noncarrier = TRUE` when the carrier table covers the analysis cohort and an omitted sample–gene pair means noncarrier. The VCF example tests one combined `DEL+DUP` mask and selects genes for follow-up by the batch q-value. Each follow-up uses the gene test's samples and covariates. It returns the alternate-allele dosage coefficient, confidence interval, carrier count, and nominal p-value. For a full single-gene fit from `run_sv_cds_test()`, call `run_sv_set_followup()` directly.

Alternatively, set `test = "burden"` for the burden test or `test = "SKAT"` for SKAT. With `keep_carrier_scores = TRUE`, `extract_sv_carrier_zscores(results, gene_id = target_gene_id)` returns scores and variant IDs for plots.

## Verification

The repository includes synthetic R tests, a two-chromosome PLINK smoke test, and WDL checks for typed file localization. See the [verification commands](docs/terra_sv_skato.md#verification). The complete workflow has yet to run on Terra.
