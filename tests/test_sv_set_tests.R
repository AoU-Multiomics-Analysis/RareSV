source("R/sv_set_tests.R")

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

test_extract_sv_by_id <- function() {
  vcf <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(vcf))
  header <- c("##fileformat=VCFv4.3", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1")
  records <- c(
    "chr20\t50\tdup.1+extra\tN\t<DUP>\t.\tPASS\tEND=125\tGT\t0/0",
    "chr20\t100\tother;dup.1+\tN\t<DUP>\t.\tPASS\tEND=120\tGT:DP\t0/1:7",
    "chr20\t150\tliteral\tCATG\tC\t.\tPASS\tSVTYPE=DEL\tGT\t0/1",
    "chr20\t200\tunknown_end\tN\t<DEL>\t.\tPASS\tSVLEN=-100\tGT\t0/0",
    "chr20\t250\tmulti\tG\tA,C\t.\tPASS\t.\tGT\t1/2"
  )
  con <- gzfile(vcf, "wt")
  writeLines(c(header, records), con)
  close(con)
  hit <- extract_sv_by_id(vcf, "dup.1+")
  stopifnot(identical(hit$record$ID, "other;dup.1+"), identical(hit$record$s1, "0/1:7"),
            identical(hit$coordinates$CHROM, "chr20"), hit$coordinates$POS == 100L,
            hit$coordinates$INFO_END == 120L, hit$coordinates$affected_start == 101L,
            hit$coordinates$affected_end == 120L)
  literal <- extract_sv_by_id(vcf, "literal")$coordinates
  stopifnot(is.na(literal$INFO_END), literal$affected_start == 151L, literal$affected_end == 153L)
  unresolved <- extract_sv_by_id(vcf, "unknown_end")$coordinates
  stopifnot(is.na(unresolved$INFO_END), is.na(unresolved$affected_end),
            extract_sv_by_id(vcf, "multi")$record$s1 == "1/2")
  expect_error(extract_sv_by_id(vcf, "dup.1"))
  expect_error(extract_sv_by_id(vcf, "."))
  con <- gzfile(vcf, "wt")
  writeLines(c(header, records, records[2]), con)
  close(con)
  expect_error(extract_sv_by_id(vcf, "dup.1+"))
}
test_extract_sv_by_id()

stopifnot(
  identical(make_sv_region("chr15", 73683944), "chr15:72683944-74683944"),
  identical(make_sv_region("chr1", 100, 1000), "chr1:1-1100")
)
expect_error(make_sv_region("chr1", NA_real_))
expect_error(make_sv_region("chr1", 10.5))

tss_locations <- tibble(
  gene_id = c("ENSG00000290825.2", "ENSG00000223972.6", "ENSG00000310526.1"),
  gene_name = c("DDX11L16", "DDX11L1", "WASH7P"),
  gene_type = c("lncRNA", "transcribed_unprocessed_pseudogene", "lncRNA"),
  chrom = rep("chr1", 3), strand = c("+", "+", "-"),
  tss_1based = c(11121L, 12010L, 30744L)
)
location <- lookup_gene_tss(tss_locations, "ENSG00000290825.2", id_col = "gene_id")
stopifnot(nrow(location) == 1L, identical(location$chrom, "chr1"), location$tss_1based == 11121)
negative <- lookup_gene_tss(tss_locations, "WASH7P", id_col = "gene_name")
stopifnot(negative$strand == "-", make_sv_region(negative$chrom, negative$tss_1based, 10) == "chr1:30734-30754")
stopifnot(nrow(lookup_gene_tss(bind_rows(tss_locations, tss_locations), "WASH7P", "gene_name")) == 1L)
expect_error(lookup_gene_tss(tss_locations, "absent"))
expect_error(lookup_gene_tss(tss_locations, "ENSG00000290825.9"))
stopifnot(lookup_gene_tss(tss_locations, "ENSG00000290825.9", ignore_gene_version = TRUE)$tss_1based == 11121)
ambiguous <- bind_rows(tss_locations, mutate(tss_locations[1, ], tss_1based = 20000L))
expect_error(lookup_gene_tss(ambiguous, "ENSG00000290825.2"))
bad_location <- tss_locations
bad_location$tss_1based[1] <- 0L
expect_error(lookup_gene_tss(bad_location, "ENSG00000290825.2"))
bad_location$tss_1based[1] <- NA_integer_
expect_error(lookup_gene_tss(bad_location, "ENSG00000290825.2"))
bad_location <- mutate(tss_locations, tss_1based = as.character(tss_1based))
bad_location$tss_1based[1] <- "bad"
expect_error(lookup_gene_tss(bad_location, "ENSG00000290825.2"))
tmp_tss <- tempfile(fileext = ".tsv")
readr::write_tsv(tss_locations, tmp_tss)
stopifnot(lookup_gene_tss(tmp_tss, "ENSG00000290825.2")$tss_1based == 11121)
unlink(tmp_tss)

gene_location <- .sv_resolve_location(
  tibble(gene_id = "ENSG00000290825.2", s1 = 1), "ENSG00000290825.2",
  tss_locations = tss_locations
)
stopifnot(gene_location$chromosome == "chr1", gene_location$tss == 11121,
          gene_location$lookup_col == "gene_id")
assay_location <- .sv_resolve_location(
  tibble(Assay = "WASH7P", s1 = 1), "WASH7P", tss_locations = tss_locations
)
stopifnot(assay_location$tss == 30744, assay_location$lookup_col == "gene_name")
symbol_pheno <- tibble(gene_name = "WASH7P", gene_id = "ENSG00000310526.1", s1 = 1, s2 = 2)
symbol_location <- .sv_resolve_location(
  symbol_pheno, "WASH7P", id_col = "gene_name", tss_locations = tss_locations
)
stopifnot(symbol_location$tss == 30744, symbol_location$lookup_col == "gene_name")
inferred_symbol <- .sv_resolve_location(
  dplyr::select(symbol_pheno, -gene_id), "WASH7P", tss_locations = tss_locations
)
stopifnot(inferred_symbol$tss == 30744, inferred_symbol$lookup_col == "gene_name")
want_symbol_pheno <- tibble(sample_id = c("s1", "s2"), phenotype = c(1, 2))
stopifnot(identical(prep_burden_phenotype(symbol_pheno, "WASH7P", "gene_name"), want_symbol_pheno))
stopifnot(identical(prep_burden_phenotype(symbol_pheno, "ENSG00000310526.1", "gene_id"), want_symbol_pheno))
stopifnot(identical(prep_burden_phenotype(dplyr::select(symbol_pheno, -gene_id), "WASH7P"), want_symbol_pheno))
override_location <- .sv_resolve_location(
  tibble(Assay = "protein_alias", s1 = 1), "protein_alias", tss_locations = tss_locations,
  tss_id = "ENSG00000290825.2", tss_id_col = "gene_id"
)
stopifnot(override_location$tss == 11121)
manual <- .sv_resolve_location(tibble(Assay = "target", s1 = 1), "target", "chr11", 93787153)
stopifnot(manual$chromosome == "chr11", manual$tss == 93787153, is.null(manual$annotation))
expect_error(.sv_resolve_location(tibble(Assay = "WASH7P"), "WASH7P", "chr1", 30744,
                                  tss_locations = tss_locations))
expect_error(.sv_resolve_location(tibble(Assay = "WASH7P"), "WASH7P"))

wide <- tibble(
  Assay = c("CD276", "IL6"), gene_id = c("g1", "g2"),
  s2 = c("1.25", "8"), s1 = c(-0.5, 3), s3 = c("NA", "4")
)
pheno <- suppressMessages(prep_burden_phenotype(
  wide, "CD276", id_col = "Assay", sample_ids = c("s1", "absent", "s3", "s2")
))
stopifnot(identical(pheno, tibble(sample_id = c("s1", "s2"), phenotype = c(-0.5, 1.25))))
stopifnot(identical(
  suppressMessages(prep_burden_phenotype(wide, "g1", id_col = "gene_id")),
  tibble(sample_id = c("s2", "s1"), phenotype = c(1.25, -0.5))
))
expect_error(prep_burden_phenotype(wide, "g1"))
expect_error(prep_burden_phenotype(bind_rows(wide, wide), "g1", id_col = "gene_id"))
expect_error(prep_burden_phenotype(tibble(Assay = "CD276", s1 = "bad"), "CD276"))

raw <- tibble(
  ID = c("sv_b", "sv_a", "sv_partial", "sv_noGT"),
  ALT = rep("<DEL>", 4), FORMAT = c("GT:DP", "DP:GT", "GT", "DP"),
  s1 = c("0/1:10", "20:1|1", "0/.", "9"),
  s2 = c("0/0:3", "8:./.", ".|.", "5")
)
want <- matrix(
  c(0L, 1L, 0L, 2L, 0L, NA_integer_, NA_integer_, NA_integer_),
  nrow = 2,
  dimnames = list(c("s2", "s1"), c("sv_b", "sv_a", "sv_partial", "sv_noGT"))
)
stopifnot(identical(parse_sv_genotypes(raw, c("s2", "s1")), want))
keep_missing <- parse_sv_genotypes(raw, c("s2", "s1"), missing_as_reference = FALSE)
stopifnot(is.na(keep_missing["s2", "sv_a"]), is.na(keep_missing["s2", "sv_partial"]))
one <- tibble(ID = "sv1", ALT = "<DEL>", FORMAT = "GT", s1 = "|0|1")
stopifnot(identical(parse_sv_genotypes(one, "s1"), matrix(1L, 1, dimnames = list("s1", "sv1"))))
bad <- raw
bad$s1[1] <- "0/2:10"
expect_error(parse_sv_genotypes(bad, "s1"))
bad$ID[2] <- bad$ID[1]
expect_error(parse_sv_genotypes(bad, "s1"))
expect_error(parse_sv_genotypes(raw, c("s1", "s1")))
expect_error(parse_sv_genotypes(raw, "absent"))
expect_error(parse_sv_genotypes(raw[0, ], "s1"))

mask_records <- tibble(
  ID = c("dup1", "inv1", "del1", "del2", "dup2", "sequence_del", "chr20-36948791-DEL-8982",
         "short_del", "DEL_in_id_only", "sequence_dup", "long_replacement"),
  REF = c(rep("A", 5), "ACGT", paste0("C", strrep("A", 8982L)), "CA", "T", "C", strrep("A", 51L)),
  ALT = c("<DUP>", "<INV>", "<DEL>", "<DEL:ME:ALU>", "<DUP:TANDEM>", "A", "C", "C", "C", "CAA", "C"),
  INFO = c("SVTYPE=DEL", ".", ".", ".", ".", "END=4;SVTYPE=DEL;FLAG", ".", ".", ".", "SVTYPE=DUP", ".")
)
stopifnot(
  identical(mask_sv_types(mask_records), mask_records),
  identical(mask_sv_types(mask_records, "DEL")$ID,
            c("del1", "del2", "sequence_del", "chr20-36948791-DEL-8982")),
  identical(mask_sv_types(mask_records, "DUP")$ID, c("dup1", "dup2", "sequence_dup")),
  identical(mask_sv_types(mask_records, "INV")$ID, "inv1"),
  identical(mask_sv_types(mask_records, c("DEL", "DUP"))$ID,
            c("dup1", "del1", "del2", "dup2", "sequence_del", "chr20-36948791-DEL-8982", "sequence_dup"))
)
expect_error(mask_sv_types(mask_records, "INS"))
expect_error(mask_sv_types(mask_records, character()))
expect_error(mask_sv_types(mask_records, c("DEL", NA_character_)))

z <- matrix(0L, 100, 4, dimnames = list(paste0("s", 1:100), c("rare", "common", "zero", "missing")))
z[1, "rare"] <- 1L
z[1:10, "common"] <- 1L
z[, "missing"] <- NA_integer_
filtered <- filter_sv_genotypes(z)
stopifnot(identical(dim(filtered$Z), c(100L, 1L)), identical(colnames(filtered$Z), "rare"))
stopifnot(identical(unname(filtered$variant_qc$keep), c(TRUE, FALSE, FALSE, FALSE)))
stopifnot(filtered$summary$all_genotypes_missing == 1L, filtered$summary$variants_retained == 1L)
none <- filter_sv_genotypes(z, max_alt_af = 0.001)
stopifnot(ncol(none$Z) == 0L)
expect_error(filter_sv_genotypes(z, max_missing = 2))

cov <- tibble(sample_id = c("s2", "s1"), age = c(NA_real_, 50))
prepared <- suppressMessages(.sv_model_data(pheno, cov, "age"))
stopifnot(identical(prepared$data$sample_id, "s1"), identical(all.vars(prepared$formula), c("phenotype", "age")))
expect_error(.sv_model_data(pheno, bind_rows(cov, cov), "age"))

test_lof_covariate <- function() {
  model <- suppressMessages(.sv_model_data(
    tibble(sample_id = paste0("s", 1:5), phenotype = c(-2, -1, 0, 1, 2)),
    tibble(sample_id = paste0("s", 1:5), age = c(50, 45, 30, 65, NA_real_)),
    "age"
  ))
  input <- tibble(
    sample_id = c("s1", "s2", "s3", "s1", "s1"),
    gene_id = c("ENSG00000123329.20", "ENSG00000123329.20", "ENSG00000123329.20",
                "ENSG00000123329.21", "ENSG00000123330.1"),
    has_lof_variant = c("true", "FALSE", NA, "1", "false")
  )
  lof <- .sv_prepare_lof_table(input)
  stopifnot(nrow(lof) == 4L, identical(lof$has_lof_variant, c(1L, 0L, NA_integer_, 0L)))
  expect_error(.sv_prepare_lof_table(bind_rows(input, mutate(input[1, ], has_lof_variant = "0"))))
  expect_error(.sv_prepare_lof_table(mutate(input, has_lof_variant = "unexpected")))

  status <- filter(lof, gene_id == "ENSG00000123329")
  known <- suppressMessages(.sv_with_lof_covariate(model, status, FALSE))
  stopifnot(
    identical(known$data$sample_id, c("s1", "s2")),
    identical(known$data$has_lof_variant, c(1L, 0L)),
    identical(attr(terms(known$formula), "term.labels"), c("age", "has_lof_variant")),
    known$lof_qc$n_lof_missing == 2L,
    known$lof_qc$n_lof_carriers == 1L,
    known$lof_qc$lof_covariate_used
  )
  sparse <- suppressMessages(.sv_with_lof_covariate(model, status, TRUE))
  stopifnot(
    identical(sparse$data$sample_id, c("s1", "s2", "s4")),
    identical(sparse$data$has_lof_variant, c(1L, 0L, 0L)),
    sparse$lof_qc$n_lof_missing == 1L
  )
  constant <- suppressMessages(.sv_with_lof_covariate(model, status[0, ], TRUE))
  stopifnot(
    identical(constant$data$sample_id, model$data$sample_id),
    identical(constant$formula, model$formula),
    constant$lof_qc$n_lof_carriers == 0L,
    !constant$lof_qc$lof_covariate_used
  )
}
test_lof_covariate()

test_lof_analysis_samples <- function() {
  test_env <- new.env()
  sys.source("R/sv_set_tests.R", envir = test_env)
  gene_id <- "ENSG00000123329"
  expression <- tibble(gene_id = gene_id, s1 = -1, s2 = 2, s3 = 0, s4 = 3, expression_only = 4)
  raw <- tibble(
    CHROM = "chr12", POS = 100L, ID = "del1", REF = "C", ALT = "<DEL>",
    QUAL = ".", FILTER = "PASS", INFO = "END=150", FORMAT = "GT",
    s1 = "0/1", s2 = "0/0", s3 = "0/1", s4 = "0/0", sv_only = "0/1"
  )
  lof <- tibble(
    sample_id = c("s1", "s3", "s4", "sv_only", "expression_only", "lof_only"),
    gene_id = paste0(gene_id, ".20"),
    has_lof_variant = c(TRUE, NA, FALSE, TRUE, TRUE, TRUE)
  )
  filtered <- NULL
  test_env$filter_sv_genotypes <- function(Z, ...) {
    stopifnot(identical(rownames(Z), c("s1", "s2", "s4")))
    filtered <<- filter_sv_genotypes(Z, ...)
    filtered
  }
  result <- tryCatch(
    suppressMessages(suppressWarnings(test_env$.sv_fit_set(
      raw, expression, gene_id, "gene_id", "DEL", "chr12:100-150", "SKATO",
      max_alt_af = 0.01, max_missing = 0.05, missing_as_reference = TRUE,
      covariates = NULL, covariate_cols = character(), n_svs_extracted = 1L,
      lof_status = .sv_prepare_lof_table(lof), lof_missing_as_noncarrier = TRUE
    ))),
    error = conditionMessage
  )
  stopifnot(
    isTRUE(all.equal(unname(filtered$variant_qc$alt_af), 1 / 6)),
    is.character(result), grepl("No SVs pass the filters", result, fixed = TRUE)
  )
}
test_lof_analysis_samples()

tmp <- tempfile(fileext = ".tsv")
writeLines(c("Assay\ts1\ts2", "CD276\t1\t2"), tmp)
stopifnot(identical(names(read_phenotype_table(tmp, "Assay")), c("Assay", "s1", "s2")))
writeLines(c("V1\tV2\tV3", "Assay\ts1\ts2", "CD276\t1\t2"), tmp)
stopifnot(identical(names(read_phenotype_table(tmp, "Assay")), c("Assay", "s1", "s2")))
unlink(tmp)

test_sv_set_batch <- function() {
  # Exercise batch behavior with deterministic outputs from the single-gene test.
  batch_env <- new.env()
  sys.source("R/sv_set_tests.R", envir = batch_env)
  batch_env$requireNamespace <- function(package, quietly) identical(package, "SKAT")
  annotation_file <- tempfile(fileext = ".tsv")
  vcf_file <- tempfile(fileext = ".vcf.gz")
  readr::write_tsv(tss_locations, annotation_file)
  file.create(vcf_file)
  on.exit(unlink(c(annotation_file, vcf_file)))
  batch_env$run_sv_set_test <- function(phenotype_df, phenotype_id, vcf_file,
                                        tss_locations, id_col, tabix_path, test,
                                        flank, sv_types) {
    stopifnot(
      identical(phenotype_df$gene_id, phenotype_id), is.data.frame(tss_locations),
      identical(sv_types, c("DEL", "DUP")), identical(flank, 100000L),
      identical(test, "SKATO")
    )
    if (phenotype_id == "no_svs") stop("No SVs pass the filters.")
    if (phenotype_id == "invalid_p") {
      return(list(summary = tibble(phenotype_id = phenotype_id, p_value = NA_real_)))
    }
    if (phenotype_id == "g1") warning("One SV remains.")
    list(summary = tibble(phenotype_id = phenotype_id, p_value = if (phenotype_id == "g1") 0.01 else 0.03))
  }
  results <- batch_env$run_sv_set_tests(
    tibble(gene_id = c("g2", "no_svs", "g1", "dup", "dup", "invalid_p"), s1 = 1:6),
    vcf_file = vcf_file, tss_locations = annotation_file, tabix_path = Sys.which("true"),
    flank = 100000L, sv_types = c("DEL", "DUP")
  )
  stopifnot(
    identical(results$phenotype_id, c("g2", "no_svs", "g1", "dup", "invalid_p")),
    identical(results$status, c("ok", "error", "ok", "error", "error")),
    identical(results$p_value, c(0.03, NA_real_, 0.01, NA_real_, NA_real_)),
    identical(results$q_value, c(0.03, NA_real_, 0.02, NA_real_, NA_real_)),
    grepl("No SVs pass", results$error_message[2]),
    grepl("multiple phenotype rows", results$error_message[4]),
    grepl("valid p-value", results$error_message[5]),
    identical(results$warnings[3], "One SV remains.")
  )
}
test_sv_set_batch()

test_targeted_sv_associations <- function() {
  test_env <- new.env()
  sys.source("R/sv_set_tests.R", envir = test_env)
  vcf_file <- tempfile()
  file.create(vcf_file)
  on.exit(unlink(vcf_file))

  samples <- paste0("s", 1:12)
  dosage <- c(rep(1L, 3), rep(0L, 9))
  phenotype <- c(2.6, 1.8, 2.4, -1.1, -1.5, -2, -.1, .3, .1, .5, .2, -.3)
  phenotype_df <- bind_cols(
    tibble(gene_id = tss_locations$gene_id[1:2]),
    as_tibble(as.list(setNames(phenotype, samples)))
  )
  gt <- rbind(ifelse(dosage == 1L, "0/1", "./."), rep("0/2", length(samples)))
  colnames(gt) <- samples
  records <- bind_cols(
    tibble(CHROM = "chr1", POS = c("11121", "11122"), ID = c("target", "unrequested"),
           REF = "N", ALT = "<DEL>", QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT"),
    as_tibble(gt)
  )
  queries <- character()
  test_env$extract_sv_region <- function(vcf_file, region, tabix_path, allow_empty) {
    queries <<- c(queries, region)
    stopifnot(allow_empty)
    if (region == "chr1:11111-11131") records else records[0, ]
  }
  targets <- tibble(
    phenotype_id = tss_locations$gene_id[c(1, 2, 1, 1)],
    variant_id = c("missing_id", "target", "target", "target")
  )
  covariates <- tibble(sample_id = rev(samples), PC1 = rep(1:4, 3))
  out <- test_env$run_targeted_sv_tests(
    targets, phenotype_df, vcf_file, tss_locations, flank = 10L,
    tabix_path = Sys.which("true"), covariates = covariates, covariate_cols = "PC1"
  )
  ref <- summary(lm(phenotype ~ rev(covariates$PC1) + dosage))$coefficients
  stopifnot(
    identical(queries, c("chr1:11111-11131", "chr1:12000-12020")),
    identical(out$variant_id, c("missing_id", "target", "target")),
    identical(out$status, c("not_in_window", "not_in_window", "ok")),
    identical(out$carriers, c(NA_integer_, NA_integer_, 3L)),
    isTRUE(all.equal(out$beta[3], unname(ref[3, 1]))),
    isTRUE(all.equal(out$p_value[3], unname(ref[3, 4])))
  )
  rare_only <- test_env$run_targeted_sv_tests(
    targets[3, ], phenotype_df, vcf_file, tss_locations, flank = 10L,
    tabix_path = Sys.which("true"), max_alt_af = 0.01
  )
  stopifnot(rare_only$status == "filtered", is.na(rare_only$p_value), rare_only$alt_af == 0.125)
}
test_targeted_sv_associations()

test_sv_set_followup <- function() {
  sample_ids <- paste0("s", seq_len(8L))
  dosage <- c(1L, 1L, 0L, 0L, NA_integer_, 0L, 0L, 0L)
  Z <- cbind(used = dosage, dropped = rev(dosage))
  rownames(Z) <- sample_ids
  data <- tibble(sample_id = sample_ids, phenotype = c(3, 2, 0, -1, .5, 0, -1, .5),
                 PC1 = rep(1:4, 2))
  result <- list(
    summary = tibble(phenotype_id = "g1", gene_id = "g1", gene_symbol = "GENE1",
                     region = "chr11:90-110", n_svs_tested = 1L),
    variants = tibble(ID = c("dropped", "used"), CHROM = "chr11", POS = c(98L, 100L), ALT = "<DEL>"),
    phenotype = data, Z_filtered = Z,
    test_result = list(test.snp.mac = c(used = 2)),
    settings = list(null_formula = "phenotype ~ PC1", max_missing = 0.25)
  )
  out <- run_sv_set_followup(result)
  data$dosage <- coalesce(dosage, mean(dosage, na.rm = TRUE))
  ref <- summary(lm(phenotype ~ PC1 + dosage, data = data))$coefficients
  stopifnot(
    identical(out$variant_id, "used"), identical(out$gene_symbol, "GENE1"),
    !"region" %in% names(out),
    identical(out$status, "ok"), identical(out$POS, 100L), identical(out$n_samples, 8L),
    isTRUE(all.equal(out$beta, unname(ref["dosage", "Estimate"]))),
    isTRUE(all.equal(out$p_value, unname(ref["dosage", "Pr(>|t|)"])))
  )
  result$test_result$test.snp.mac <- c(other = 2)
  expect_error(run_sv_set_followup(result))
  result$test_result$test.snp.mac <- c(2)
  expect_error(run_sv_set_followup(result))
  result$Z_filtered[1L, "dropped"] <- NA_integer_
  stopifnot(identical(run_sv_set_followup(result)$variant_id, "used"))
  result$Z_filtered <- result$Z_filtered[, "used", drop = FALSE]
  stopifnot(identical(run_sv_set_followup(result)$variant_id, "used"))
}
test_sv_set_followup()

test_restore_sv_followup_ids <- function() {
  env <- new.env()
  sys.source("R/sv_set_tests.R", envir = env)
  records <- tibble(
    CHROM = "chr12", POS = c("124912413", "124912413", "124912300"),
    ID = c("chr12-124912413-DEL-60", "chr12-124912413-DEL-72", "spanning"),
    REF = c(str_c("A", str_dup("C", 60L)), str_c("A", str_dup("G", 72L)), "N"),
    ALT = c("A", "A", "<DEL>"), INFO = c(".", ".", "END=124912500")
  )
  queries <- character()
  env$extract_sv_region <- function(vcf_file, region, tabix_path, allow_empty, include_samples) {
    stopifnot(!include_samples, allow_empty)
    queries <<- c(queries, region)
    records
  }
  old <- tibble(gene_symbol = "UBC", variant_id = c("12_206463", "old", records$ID[2], "elsewhere"),
                region = "chr12:124912400-124912500",
                CHROM = "chr12", POS = c(124912413L, 124912413L, 124912413L, 124912414L),
                ALT = "A", REF = c(records$REF[1], NA_character_, NA_character_, NA_character_),
                beta = c(1.51, 1.1, .8, 1))
  out <- env$restore_sv_followup_ids(old, "original.vcf.gz")
  stopifnot(!any(c("region", "previous_variant_id") %in% names(out)), identical(out$beta, old$beta),
            identical(out$variant_id, c(records$ID[1], NA_character_, records$ID[2], NA_character_)),
            identical(out$id_match_status, c("matched", "ambiguous", "matched", "not_found")),
            length(unique(queries)) == 2L, length(queries) == 2L)
  records <- records[-2, ]
  stopifnot(identical(env$restore_sv_followup_ids(dplyr::select(old[1, ], -REF, -variant_id, -region), "original.vcf.gz")$variant_id,
                      records$ID[1]))
}
test_restore_sv_followup_ids()

test_sv_carrier_zscores <- function() {
  sample_ids <- paste0("s", 1:5)
  Z <- cbind(
    del1 = c(1L, 0L, NA_integer_, 0L, 0L),
    del2 = c(2L, 1L, 0L, 0L, 0L),
    dup1 = c(1L, 0L, 0L, 1L, NA_integer_),
    dropped = c(0L, 0L, 1L, 0L, 1L),
    inv1 = c(0L, 0L, 1L, 0L, 0L)
  )
  rownames(Z) <- sample_ids
  result <- list(
    summary = tibble(phenotype_id = "g1", gene_id = "g1", gene_symbol = "GENE1", n_svs_tested = 4L),
    phenotype = tibble(sample_id = sample_ids, phenotype = c(-2, -1, 0, 1, 2)),
    Z_filtered = Z,
    variants = tibble(
      ID = c("dup1", "dropped", "del2", "del1", "inv1"),
      ALT = c("AA", "<DEL>", "<DEL:ME>", "<DEL>", "<INV>"),
      INFO = c("SVTYPE=DUP", ".", ".", ".", ".")
    ),
    test_result = list(test.snp.mac = c(del1 = 1, del2 = 3, dup1 = 2, inv1 = 1))
  )
  expected <- tibble(
    phenotype_id = "g1", gene_id = "g1", gene_symbol = "GENE1",
    sample_id = c("s1", "s2", "s1", "s4"), sv_type = c("DEL", "DEL", "DUP", "DUP"),
    z_score = c(-2, -1, -2, 1), n_carried_svs = c(2L, 1L, 1L, 1L),
    variant_ids = c("del1,del2", "del2", "dup1", "dup1")
  )
  stopifnot(identical(extract_sv_carrier_zscores(result), expected))

  result$summary$n_svs_tested <- 1L
  result$test_result$test.snp.mac <- c(inv1 = 1)
  stopifnot(identical(extract_sv_carrier_zscores(result), expected[0, ]))
  result$phenotype <- result$phenotype[5:1, ]
  expect_error(extract_sv_carrier_zscores(result))
}
test_sv_carrier_zscores()

test_sv_sheet_associations <- function() {
  test_env <- new.env()
  sys.source("R/sv_set_tests.R", envir = test_env)
  vcf_file <- tempfile()
  file.create(vcf_file)
  on.exit(unlink(vcf_file))

  samples <- paste0("s", 1:8)
  d <- c(1L, 1L, rep(0L, 6))
  y <- c(2, 3, 0, -1, .5, 0, -1, .5)
  pheno <- bind_cols(tibble(gene_id = c("g1", "g2"), gene_symbol = c("GENE1", "GENE2")),
                    as_tibble(as.list(setNames(y, samples))))
  gt <- rbind(ifelse(d == 1L, "0/1", "./."), rep("0/2", 8), rep("0/2", 8),
              ifelse(d == 1L, "0/1", "./."))
  colnames(gt) <- samples
  records <- bind_cols(
    tibble(CHROM = "chr11", POS = c("100", "100", "99", "200"),
           ID = c("sv1", "sv2", "spanning_sv", "."),
           REF = "N", ALT = "<DEL>", QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT"),
    as_tibble(gt)
  )
  queries <- character()
  test_env$extract_sv_region <- function(vcf_file, region, tabix_path, allow_empty) {
    queries <<- c(queries, region)
    stopifnot(allow_empty)
    if (region == "chr11:100-100") return(records[1:3, ])
    if (region == "chr11:200-200") return(records[4, ])
    records[0, ]
  }
  sheet <- tibble(
    Ensembl = c("g2", "g1", "g1", "g1", "g1"), `SV chromosome` = "chr11",
    `VCF position` = c("100", "200", "100", "300", "100"),
    `VCF ID` = c("sv1", NA, NA, "sv1", "spanning_sv")
  )
  out <- test_env$run_sv_sheet_tests(
    sheet, pheno, vcf_file, gene_col = "Ensembl", chrom_col = "SV chromosome",
    pos_col = "VCF position", variant_id_col = "VCF ID", tabix_path = Sys.which("true")
  )
  ref <- summary(lm(y ~ d))$coefficients
  stopifnot(
    identical(sort(queries), c("chr11:100-100", "chr11:200-200", "chr11:300-300")),
    identical(out$sheet_row, 1:5),
    identical(out$status, c("ok", "ok", "ambiguous", "not_found", "not_found")),
    identical(out$phenotype_id, sheet$Ensembl),
    identical(out$gene_id, sheet$Ensembl),
    identical(out$gene_symbol, c("GENE2", rep("GENE1", 4))),
    identical(out$variant_id[1:2], c("sv1", NA_character_)),
    grepl("sv1", out$error_message[3]), grepl("sv2", out$error_message[3]),
    !grepl("spanning_sv", out$error_message[3]),
    isTRUE(all.equal(out$beta[1:2], rep(unname(ref[2, 1]), 2))),
    isTRUE(all.equal(out$p_value[1:2], rep(unname(ref[2, 4]), 2)))
  )
  coordinate_only <- test_env$run_sv_sheet_tests(
    dplyr::select(sheet[2, ], -`VCF ID`), pheno, vcf_file, gene_col = "Ensembl",
    chrom_col = "SV chromosome", pos_col = "VCF position", tabix_path = Sys.which("true")
  )
  stopifnot(coordinate_only$status == "ok", is.na(coordinate_only$variant_id))
}
test_sv_sheet_associations()

test_sv_cds_sets <- function() {
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    message("Skipping the CDS fixture because IRanges is unavailable.")
    return(invisible(NULL))
  }
  test_env <- new.env()
  sys.source("R/sv_set_tests.R", envir = test_env)
  gtf <- tempfile(fileext = ".gtf")
  on.exit(unlink(gtf))
  writeLines(c(
    'chr20\ttest\tCDS\t100\t110\t.\t-\t0\tgene_id "g1.2"; transcript_id "t1";',
    'chr20\ttest\tCDS\t105\t112\t.\t-\t0\tgene_id "g1.2"; transcript_id "t2";',
    'chr20\ttest\tCDS\t200\t210\t.\t-\t0\tgene_id "g1.2"; transcript_id "t1";',
    'chr20\ttest\texon\t140\t180\t.\t-\t.\tgene_id "g1.2"; transcript_id "t1";',
    'chr20\ttest\tCDS\t140\t180\t.\t+\t0\tgene_id "g2"; transcript_id "t3";'
  ), gtf)
  raw <- tibble(
    CHROM = "chr20", POS = c(89L, 112L, 150L, 199L, 190L, 200L),
    ID = c("upstream_del", "anchor_only", "intronic", "dup", "literal_del", "no_span"),
    REF = c(rep("C", 4), paste0("C", strrep("A", 55L)), "C"),
    ALT = c("<DEL>", "<DEL>", "<DEL>", "<DUP>", "C", "<DEL>"),
    INFO = c("END=112", "END=120", "END=170", "END=201", ".", "SVLEN=-50"),
    s1 = "0/1", s2 = "0/0"
  )
  test_env$extract_sv_region <- function(vcf_file, region, tabix_path, allow_empty) {
    stopifnot(region == "chr20:100-210", allow_empty)
    raw
  }
  cds <- test_env$read_cds_gtf(gtf)
  out <- test_env$extract_cds_svs("unused.vcf.gz", cds, "g1", ignore_gene_version = TRUE)
  stopifnot(
    identical(out$region, "chr20:100-210"),
    identical(out$cds_intervals$start, c(100L, 200L)),
    identical(out$cds_intervals$end, c(112L, 210L)),
    identical(out$variants$ID, c("upstream_del", "dup", "literal_del")),
    identical(out$overlap$cds_overlap_bp, c(13, 0, 0, 2, 11, NA_real_)),
    identical(out$overlap$affected_start, c(90L, 113L, 151L, 200L, 191L, NA_integer_))
  )
  expect_error(test_env$extract_cds_svs("unused.vcf.gz", cds, "g1"))
  expect_error(test_env$extract_cds_svs("unused.vcf.gz", cds, "gene_with_no_cds"))

  test_env$requireNamespace <- function(package, quietly) {
    if (package == "SKAT") TRUE else base::requireNamespace(package, quietly = quietly)
  }
  test_env$.sv_fit_set <- function(raw_sv, phenotype_df, phenotype_id, id_col, sv_types,
                                   region, test, max_alt_af, max_missing, missing_as_reference,
                                   covariates, covariate_cols, n_svs_extracted, ...) {
    stopifnot(identical(raw_sv$ID, c("upstream_del", "dup", "literal_del")),
              n_svs_extracted == 6L, phenotype_id == "g1", test == "burden",
              identical(sv_types, c("DEL", "DUP")))
    list(summary = tibble(region = region, p_value = 0.01), settings = list(test = test))
  }
  tested <- test_env$run_sv_cds_test(tibble(gene_id = "g1"), "g1", "unused.vcf.gz", cds,
                                     ignore_gene_version = TRUE, test = "burden")
  stopifnot(tested$summary$n_svs_overlapping_cds == 3L,
            tested$summary$n_svs_unresolved_span == 1L,
            tested$settings$coordinate_source == "gtf_cds")
  test_env$extract_sv_region <- function(vcf_file, region, tabix_path) {
    stopifnot(region == "chr20:100-210")
    raw
  }
  test_env$.sv_fit_set <- function(raw_sv, ..., n_svs_extracted) {
    stopifnot(identical(raw_sv$ID, raw$ID), n_svs_extracted == nrow(raw))
    list(summary = tibble(p_value = 0.01), settings = list(test = "SKATO"))
  }
  legacy <- test_env$run_sv_set_test(tibble(gene_id = "g1"), "g1", "unused.vcf.gz",
                                    "chr20", 155L, flank = 55L)
  stopifnot("tss_location" %in% names(legacy), is.null(legacy$tss_location),
            legacy$settings$tss == 155L, legacy$settings$coordinate_source == "manual")
}
test_sv_cds_sets()

cat("PASS: TSS and CDS selection, SV windows, phenotype and genotype preparation, ALT masks, variant filters, batch summaries, and targeted SV tests.\n")
