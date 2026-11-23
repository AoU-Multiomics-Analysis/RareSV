source("R/sv_set_tests.R")

g1 <- "ENSG00000101347"
g2 <- "ENSG00000100002"
reference <- tibble(
  gene_id = c(paste0(g1, ".4"), paste0(g2, ".1"), "ENSG00000100003.2", "ENSG00000100003.2_PAR_Y"),
  gene_name = c("SAMHD1", "GENE2", "PAR", "PAR"),
  chrom = c("chr20", "chr2", "chrX", "chrY"), tss_1based = 100L, strand = "+"
)
cds <- tibble(gene_id = paste0(g1, ".7"), chrom = "chr20", start = 105L, end = 115L)
sheet <- tibble(gene_id = c(paste0(g1, ".8"), NA, "PAR", g1, "unknown"),
                gene_name = c("SAMHD1", "GENE2", NA, "GENE2", NA), pos = 1:5)
expr <- tibble(target = c("SAMHD1", paste0(g2, ".4"), g2), s1 = 1:3)
protein <- tibble(Assay = c("SAMHD1", "SAMHD1_extra"), s1 = c(1, 2))
h <- harmonize_sv_inputs(
  reference, sv_sheet = sheet, expression = expr, proteomics = protein,
  cds_gtf = cds, expression_col = "target",
  assay_map = tibble(Assay = "SAMHD1_extra", gene_id = g1)
)
stopifnot(
  identical(h$sv_sheet$gene_id, c(g1, g2)), identical(h$sv_sheet$pos, 1:2),
  identical(h$sv_sheet$source_sheet_row, 1:2),
  identical(h$expression$gene_id, g1), identical(h$expression$s1, 1L),
  identical(h$proteomics$Assay, protein$Assay), identical(h$proteomics$gene_id, rep(g1, 2)),
  identical(h$cds_gtf$gene_id, g1), identical(h$cds_gtf$gtf_gene_id, paste0(g1, ".7")),
  "ENSG00000100003_PAR_Y" %in% h$tss_locations$gene_id,
  identical(h$mapping_qc$status[h$mapping_qc$dataset == "sv_sheet"],
            c("ok", "ok", "ambiguous", "conflict", "unmapped")),
  identical(h$mapping_qc$status[h$mapping_qc$dataset == "expression"],
            c("ok", "duplicate_phenotype", "duplicate_phenotype"))
)
annotated <- annotate_sv_results(tibble(phenotype_id = protein$Assay, p_value = c(.02, .03)),
                                 h$proteomics, "Assay")
stopifnot(identical(annotated$gene_id, rep(g1, 2)), identical(annotated$gene_symbol, rep("SAMHD1", 2)))
resolved <- .sv_resolve_location(h$proteomics, "SAMHD1_extra", id_col = "Assay",
                                 tss_locations = h$tss_locations)
stopifnot(resolved$chromosome == "chr20", resolved$lookup_col == "gene_id")
explicit_lookup <- .sv_resolve_location(h$proteomics, "SAMHD1_extra", id_col = "Assay",
                                        tss_locations = h$tss_locations, tss_id_col = "gene_id")
stopifnot(explicit_lookup$query_id == g1)
stopifnot(identical(prep_burden_phenotype(h$expression, g1)$sample_id, "s1"))
labelled_summary <- annotate_sv_results(list(summary = tibble(phenotype_id = "SAMHD1_extra")),
                                        h$proteomics, id_col = "Assay")$summary
stopifnot(labelled_summary$gene_id == g1, labelled_summary$gene_symbol == "SAMHD1")

test_env <- new.env()
sys.source("R/sv_set_tests.R", test_env)
test_env$requireNamespace <- function(package, quietly) identical(package, "SKAT")
test_env$extract_cds_svs <- function(vcf_file, gtf, gene_id, tabix_path, ignore_gene_version) {
  stopifnot(identical(gene_id, g1), identical(gtf$gene_id, g1))
  list(variants = tibble(ID = "sv1"), overlap = tibble(span_status = "ok", sv_type = "DEL"),
       region = "chr20:105-115", cds_intervals = tibble(start = 105L, end = 115L), gene_id = gene_id)
}
test_env$.sv_fit_set <- function(raw_sv, phenotype_df, phenotype_id, id_col, ...) {
  stopifnot(phenotype_id == "SAMHD1_extra", id_col == "Assay")
  list(summary = annotate_sv_results(tibble(phenotype_id = phenotype_id, p_value = 0.01),
                                     phenotype_df, id_col), settings = list())
}
protein_cds <- test_env$run_sv_cds_test(h$proteomics, "SAMHD1_extra", "unused", h$cds_gtf, id_col = "Assay")
stopifnot(protein_cds$settings$cds_gene_id == g1, protein_cds$summary$gene_symbol == "SAMHD1")

test_cds_batch <- function() {
  batch_env <- new.env()
  sys.source("R/sv_set_tests.R", envir = batch_env)
  batch_env$requireNamespace <- function(package, quietly) {
    if (package == "SKAT") TRUE else base::requireNamespace(package, quietly = quietly)
  }
  gtf_file <- tempfile(fileext = ".gtf")
  vcf_file <- tempfile(fileext = ".vcf.gz")
  writeLines(c(
    sprintf('chr20\ttest\tCDS\t101\t109\t.\t+\t0\tgene_id "%s.2";', g1),
    sprintf('chr2\ttest\tCDS\t101\t109\t.\t+\t0\tgene_id "%s.6";', g2)
  ), gtf_file)
  file.create(vcf_file)
  on.exit(unlink(c(gtf_file, vcf_file)))

  g3 <- "ENSG00000100003"
  absent <- "ENSG00000100003_PAR_Y"
  expression <- tibble(gene_id = c(g1, g2, g3), s1 = 1:3)
  h <- batch_env$harmonize_sv_inputs(reference, expression = expression)
  read_count <- 0L
  regions <- character()
  batch_env$read_cds_gtf <- function(path) {
    read_count <<- read_count + 1L
    read_cds_gtf(path)
  }
  batch_env$extract_sv_region <- function(vcf_file, region, tabix_path, allow_empty) {
    stopifnot(region %in% c("chr2:101-109", "chr20:101-109"), allow_empty)
    regions <<- c(regions, region)
    tibble(CHROM = str_remove(region, ":.*"), POS = 100L, ID = "coding_sv", REF = "C",
           ALT = "<DEL>", INFO = "END=105")
  }
  saved_samples <- tibble(sample_id = paste0("s", 1:8), phenotype = c(3, 2, 0, -1, .5, 0, -1, .5))
  saved_genotypes <- cbind(coding_sv = c(1L, 1L, rep(0L, 6L)))
  rownames(saved_genotypes) <- saved_samples$sample_id
  batch_env$.sv_fit_set <- function(raw_sv, phenotype_df, phenotype_id, id_col, sv_types,
                                    region, test, max_alt_af, ...) {
    stopifnot(identical(raw_sv$ID, "coding_sv"), identical(phenotype_df$gene_id, phenotype_id),
              identical(sv_types, "DEL"), identical(test, "SKATO"), identical(max_alt_af, 0.02))
    list(
      summary = tibble(phenotype_id = phenotype_id, region = region, test = test, n_svs_tested = 1L,
                       p_value = if (phenotype_id == g1) 0.01 else 0.44),
      phenotype = saved_samples, Z_filtered = saved_genotypes, variants = raw_sv,
      settings = list(null_formula = "phenotype ~ 1"), test_result = list(test.snp.mac = c(coding_sv = 2))
    )
  }
  out <- batch_env$run_sv_cds_tests(
    h$expression, c(g2, g3, paste0(g1, ".7"), g1, absent), vcf_file, gtf_file,
    tabix_path = Sys.which("true"), ignore_gene_version = TRUE, sv_types = "DEL", max_alt_af = 0.02
  )
  storey <- qvalue::qvalue(c(0.44, 0.01), lambda = seq(0.05, 0.40, 0.05), lfdr.out = FALSE)$qvalues
  stopifnot(
    identical(read_count, 1L), identical(out$gene_id, c(g2, g3, g1, absent)),
    identical(out$gene_symbol, c("GENE2", "PAR", "SAMHD1", NA_character_)),
    identical(out$status, c("ok", "error", "ok", "error")),
    identical(out$n_svs_overlapping_cds, c(1L, NA_integer_, 1L, NA_integer_)),
    isTRUE(all.equal(out$q_value, c(storey[1L], NA_real_, storey[2L], NA_real_))),
    grepl("No CDS rows match", out$error_message[2L]),
    grepl("No phenotype row", out$error_message[4L])
  )

  h <- harmonize_sv_inputs(reference, expression = expression, cds_gtf = read_cds_gtf(gtf_file))
  out <- suppressWarnings(
    batch_env$run_sv_cds_tests(
      h$expression[-2L, ], c(g2, g1), vcf_file, h$cds_gtf,
      tabix_path = Sys.which("true"), sv_types = "DEL", max_alt_af = 0.02, follow_up = TRUE
    )
  )
  stopifnot(
    identical(read_count, 1L), identical(out$summary$gene_symbol, c("GENE2", "SAMHD1")),
    identical(out$summary$status, c("error", "ok")), all(is.na(out$summary$q_value)),
    nrow(out$single_variant) == 0L, !"region" %in% names(out$single_variant)
  )

  followup_ids <- character()
  batch_env$run_sv_set_followup <- function(result) {
    id <- result$summary$phenotype_id
    followup_ids <<- c(followup_ids, id)
    if (id == g2) stop("Missing SKAT variant IDs.")
    run_sv_set_followup(result)
  }
  regions <- character()
  batch <- batch_env$run_sv_cds_tests(
    h$expression, c(g2, g3, g1), vcf_file, h$cds_gtf,
    tabix_path = Sys.which("true"), sv_types = "DEL", max_alt_af = 0.02,
    follow_up = TRUE, q_cutoff = 0.018
  )
  stopifnot(
    identical(batch$summary$status, c("ok", "error", "ok")),
    isTRUE(all.equal(batch$summary$q_value, c(storey[1L], NA_real_, storey[2L]))),
    identical(regions, c("chr2:101-109", "chr20:101-109")), identical(followup_ids, g1),
    identical(batch$single_variant$gene_id, g1),
    identical(batch$single_variant$gene_symbol, "SAMHD1"),
    identical(batch$single_variant$variant_id, "coding_sv"),
    identical(batch$single_variant$status, "ok"),
    is.finite(batch$single_variant$p_value)
  )

  original_vcf_file <- tempfile(fileext = ".vcf.gz")
  file.create(original_vcf_file)
  on.exit(unlink(original_vcf_file), add = TRUE)
  extract_analysis_region <- batch_env$extract_sv_region
  batch_env$extract_sv_region <- function(vcf_file, region, tabix_path, allow_empty, include_samples = TRUE) {
    if (!identical(vcf_file, original_vcf_file)) {
      return(extract_analysis_region(vcf_file, region, tabix_path, allow_empty))
    }
    stopifnot(region == "chr20:100-100", !include_samples)
    tibble(CHROM = "chr20", POS = 100L, ID = "source-del;source-alias", REF = "C", ALT = "<DEL>", INFO = "END=105")
  }
  original_batch <- batch_env$run_sv_cds_tests(
    h$expression, c(g2, g3, g1), vcf_file, h$cds_gtf,
    tabix_path = Sys.which("true"), sv_types = "DEL", max_alt_af = 0.02,
    follow_up = TRUE, q_cutoff = 0.018, original_vcf_file = original_vcf_file
  )
  stopifnot(identical(original_batch$single_variant$variant_id, "source-del;source-alias"),
            !any(c("region", "previous_variant_id") %in% names(original_batch$single_variant)),
            identical(original_batch$single_variant$beta, batch$single_variant$beta))

  carrier_batch <- batch_env$run_sv_cds_tests(
    h$expression, c(g2, g3, g1), vcf_file, h$cds_gtf, tabix_path = Sys.which("true"),
    sv_types = "DEL", max_alt_af = 0.02, keep_carrier_scores = TRUE
  )
  batch_file <- tempfile(fileext = ".rds")
  on.exit(unlink(batch_file), add = TRUE)
  saveRDS(carrier_batch, batch_file)
  carriers <- batch_env$extract_sv_carrier_zscores(batch_file, gene_id = paste0(g1, ".2"))
  stopifnot(
    identical(sort(names(carrier_batch)), c("carrier_scores", "summary")),
    identical(carrier_batch$summary, batch$summary),
    identical(carriers$gene_id, rep(g1, 2L)),
    identical(carriers$gene_symbol, rep("SAMHD1", 2L)),
    identical(carriers$sample_id, c("s1", "s2")),
    identical(carriers$sv_type, c("DEL", "DEL")),
    identical(carriers$z_score, c(3, 2)),
    identical(carriers$n_carried_svs, c(1L, 1L)),
    identical(carriers$variant_ids, c("coding_sv", "coding_sv"))
  )
  missing_scores <- tryCatch(
    batch_env$extract_sv_carrier_zscores(batch_file, gene_id = g3), error = conditionMessage
  )
  stopifnot(is.character(missing_scores), grepl(g3, missing_scores, fixed = TRUE))
  saveRDS(batch, batch_file)
  old_batch <- tryCatch(
    batch_env$extract_sv_carrier_zscores(batch_file, gene_id = g1), error = conditionMessage
  )
  stopifnot(is.character(old_batch), grepl("keep_carrier_scores", old_batch, fixed = TRUE))

  both <- batch_env$run_sv_cds_tests(
    h$expression, c(g2, g1), vcf_file, h$cds_gtf, tabix_path = Sys.which("true"),
    sv_types = "DEL", max_alt_af = 0.02, follow_up = TRUE, q_cutoff = 1, keep_carrier_scores = TRUE
  )
  stopifnot(
    identical(both$summary$status, c("ok", "ok")),
    identical(both$single_variant$gene_id, c(g2, g1)),
    identical(both$single_variant$gene_symbol, c("GENE2", "SAMHD1")),
    identical(both$single_variant$variant_id, c(NA_character_, "coding_sv")),
    identical(both$single_variant$status, c("error", "ok")),
    identical(batch_env$extract_sv_carrier_zscores(both, gene_id = g1), carriers)
  )

  lof <- tibble(
    sample_id = c("s1", "s2", "s2"),
    gene_id = paste0(c(g1, g1, g2), ".20"),
    has_lof_variant = c("true", "false", "true")
  )
  fit <- batch_env$.sv_fit_set
  lof_by_gene <- list()
  batch_env$.sv_fit_set <- function(raw_sv, phenotype_df, phenotype_id, ...,
                                    lof_status, lof_missing_as_noncarrier) {
    stopifnot(isTRUE(lof_missing_as_noncarrier))
    lof_by_gene[[phenotype_id]] <<- lof_status %>% dplyr::select(sample_id, has_lof_variant)
    fit(raw_sv, phenotype_df, phenotype_id, ...)
  }
  adjusted <- batch_env$run_sv_cds_tests(
    h$expression, c(g2, g1), vcf_file, h$cds_gtf, tabix_path = Sys.which("true"),
    sv_types = "DEL", max_alt_af = 0.02, lof_table = lof, lof_missing_as_noncarrier = TRUE
  )
  stopifnot(
    identical(adjusted$status, c("ok", "ok")),
    identical(lof_by_gene[[g1]], tibble(sample_id = c("s1", "s2"), has_lof_variant = c(1L, 0L))),
    identical(lof_by_gene[[g2]], tibble(sample_id = "s2", has_lof_variant = 1L))
  )
}
test_cds_batch()

cat("PASS: gene harmonization, assay identity, CDS batch order, failure rows, and adjusted p-values.\n")
