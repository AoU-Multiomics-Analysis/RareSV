# Source this file, then call run_sv_set_test().
# Inputs use a wide phenotype table and an indexed, bgzip-compressed SV VCF.
# Genotypes must be diploid and biallelic. Coordinates must use the same build.
library(tidyverse)

# Accept an annotation table or its TSV path. Keep gene versions unless requested.
# tss_1based already specifies the TSS for either strand; do not reverse it.
lookup_gene_tss <- function(tss_locations, target_id, id_col = c("gene_id", "gene_name"),
                            ignore_gene_version = FALSE) {
  id_col <- match.arg(id_col)
  stopifnot(is.logical(ignore_gene_version), length(ignore_gene_version) == 1L,
            !is.na(ignore_gene_version))
  if (length(target_id) != 1L || is.na(target_id) || !nzchar(as.character(target_id))) {
    stop("Supply one non-missing target identifier for the TSS lookup.")
  }
  if (is.character(tss_locations) && length(tss_locations) == 1L) {
    tss_locations <- data.table::fread(tss_locations, colClasses = "character", check.names = FALSE)
  }
  annotation <- as_tibble(tss_locations)
  required <- c("gene_id", "gene_name", "chrom", "strand", "tss_1based")
  if (!all(required %in% names(annotation))) {
    stop("Missing TSS columns: ", paste(setdiff(required, names(annotation)), collapse = ", "))
  }
  query_id <- as.character(target_id)
  keys <- as.character(annotation[[id_col]])
  if (id_col == "gene_id" && ignore_gene_version) {
    keys <- str_remove(keys, "\\.[0-9]+$")
    query_id <- str_remove(query_id, "\\.[0-9]+$")
  }
  location <- annotation[which(!is.na(keys) & keys == query_id), , drop = FALSE] %>% distinct()
  if (!nrow(location)) {
    stop("No TSS row matches ", id_col, " = ", target_id, ". Check the identifier and gene version.")
  }
  if (nrow(location) != 1L) {
    stop("Multiple TSS rows match ", target_id, ". Filter the annotation to one row for this target.")
  }
  chromosome <- as.character(location$chrom[[1]])
  strand_value <- as.character(location$strand[[1]])
  position <- suppressWarnings(as.numeric(as.character(location$tss_1based[[1]])))
  if (is.na(chromosome) || !nzchar(chromosome) || str_detect(chromosome, "[:\\s]")) {
    stop("Invalid chrom value for ", target_id, ".")
  }
  if (is.na(strand_value) || !strand_value %in% c("+", "-")) {
    stop("Invalid strand value for ", target_id, ".")
  }
  if (!is.finite(position) || position < 1 || position != floor(position)) {
    stop("Invalid tss_1based value for ", target_id, ". Expected a positive integer.")
  }
  location %>% mutate(chrom = .env$chromosome, strand = .env$strand_value,
                       tss_1based = .env$position)
}

.sv_harmonized_gene_id <- function(phenotype_df, phenotype_id, id_col) {
  if (is.null(id_col)) return(NULL)
  if (!all(c("gene_id", "gene_symbol", id_col) %in% names(phenotype_df))) return(NULL)
  ids <- as.character(phenotype_df$gene_id[as.character(phenotype_df[[id_col]]) %in% phenotype_id])
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) > 1L) stop("The phenotype maps to multiple Ensembl gene IDs.")
  if (length(ids)) ids[[1L]] else NULL
}

# Add labels to a result table or to the summary of a single set test.
annotate_sv_results <- function(result, phenotype_df, id_col = NULL) {
  if (!is.data.frame(result) && is.list(result) && is.data.frame(result$summary)) {
    result$summary <- annotate_sv_results(result$summary, phenotype_df, id_col)
    return(result)
  }
  if (!all(c("gene_id", "gene_symbol") %in% names(phenotype_df)) ||
      !"phenotype_id" %in% names(result)) return(result)
  if (is.null(id_col)) {
    ids <- names(phenotype_df)[str_to_lower(names(phenotype_df)) %in% c("assay", "gene_id", "gene_name")]
    if (length(ids) != 1L) stop("Specify id_col to annotate the result.")
    id_col <- ids[[1L]]
  }
  if (!id_col %in% names(phenotype_df)) stop("id_col must name an existing phenotype column.")
  labels <- tibble(
    phenotype_id = as.character(phenotype_df[[id_col]]),
    gene_id = as.character(phenotype_df$gene_id), gene_symbol = as.character(phenotype_df$gene_symbol)
  ) %>% distinct()
  if (anyDuplicated(labels$phenotype_id)) stop("A phenotype ID has conflicting gene labels.")
  result %>% dplyr::select(-any_of(c("gene_id", "gene_symbol"))) %>%
    left_join(labels, by = "phenotype_id") %>% relocate(gene_id, gene_symbol, .after = phenotype_id)
}

.sv_resolve_location <- function(phenotype_df, phenotype_id, chromosome = NULL, tss = NULL,
                                 id_col = NULL, tss_locations = NULL, tss_id = NULL,
                                 tss_id_col = NULL, ignore_gene_version = FALSE) {
  if (is.null(tss_locations)) {
    if (!is.null(tss_id) || !is.null(tss_id_col)) stop("TSS lookup options require tss_locations.")
    if (is.null(chromosome) || is.null(tss)) {
      stop("Supply tss_locations, or supply both chromosome and tss.")
    }
    make_sv_region(chromosome, tss, flank = 0)
    return(list(chromosome = chromosome, tss = tss, annotation = NULL,
                lookup_col = NULL, query_id = NULL))
  }
  if (!is.null(chromosome) || !is.null(tss)) {
    stop("Supply either tss_locations or manual chromosome and tss coordinates.")
  }
  if (is.null(id_col)) {
    candidates <- names(phenotype_df)[str_to_lower(names(phenotype_df)) %in% c("assay", "gene_id", "gene_name")]
    if (length(candidates) != 1L) stop("Specify id_col when the phenotype identifier column is ambiguous.")
    id_col <- candidates
  }
  if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(phenotype_df)) {
    stop("id_col must name one existing phenotype column.")
  }
  if (is.null(tss_id) && (is.null(tss_id_col) || identical(tss_id_col, "gene_id"))) {
    mapped_id <- .sv_harmonized_gene_id(phenotype_df, phenotype_id, id_col)
    if (!is.null(mapped_id)) {
      tss_id <- mapped_id
      tss_id_col <- "gene_id"
    }
  }
  if (is.null(tss_id_col)) {
    tss_id_col <- switch(str_to_lower(id_col), gene_id = "gene_id", gene_name = "gene_name",
                         assay = "gene_name", NULL)
    if (is.null(tss_id_col)) stop("Specify tss_id_col as gene_id or gene_name for this phenotype column.")
  }
  query_id <- if (is.null(tss_id)) phenotype_id else tss_id
  annotation <- lookup_gene_tss(tss_locations, query_id, tss_id_col, ignore_gene_version)
  list(chromosome = annotation$chrom[[1]], tss = annotation$tss_1based[[1]],
       annotation = annotation, lookup_col = tss_id_col, query_id = as.character(query_id))
}

make_sv_region <- function(chromosome, tss, flank = 1000000L) {
  stopifnot(
    is.character(chromosome), length(chromosome) == 1L,
    !is.na(chromosome), nzchar(chromosome),
    is.numeric(tss), length(tss) == 1L, is.finite(tss), tss >= 1, tss == floor(tss),
    is.numeric(flank), length(flank) == 1L, is.finite(flank), flank >= 0,
    flank == floor(flank)
  )
  str_c(chromosome, ":", format(max(1, tss - flank), scientific = FALSE, trim = TRUE),
        "-", format(tss + flank, scientific = FALSE, trim = TRUE))
}

read_phenotype_table <- function(file, id_col) {
  data <- data.table::fread(
    file, header = TRUE, colClasses = "character", check.names = FALSE
  ) %>% as_tibble(.name_repair = "minimal")

  if (!id_col %in% names(data)) {
    header <- as.character(unlist(slice_head(data, n = 1), use.names = FALSE))
    if (!id_col %in% header) stop("Cannot find phenotype identifier column: ", id_col)
    data <- slice(data, -1)
    names(data) <- header
  }
  as_tibble(data, .name_repair = "check_unique")
}

prep_burden_phenotype <- function(phenotype_df, phenotype_id,
                                  id_col = NULL, sample_ids = NULL) {
  phenotype_df <- as_tibble(phenotype_df)
  id_columns <- names(phenotype_df)[
    str_to_lower(names(phenotype_df)) %in% c("assay", "gene_id", "gene_name")
  ]
  if (is.null(id_col)) {
    if (length(id_columns) != 1L) stop("Specify id_col when the identifier column is ambiguous.")
    id_col <- id_columns
  }
  if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(phenotype_df)) {
    stop("id_col must name one existing column.")
  }
  if (length(phenotype_id) != 1L || is.na(phenotype_id)) {
    stop("Supply one non-missing phenotype_id.")
  }

  selected <- phenotype_df %>%
    filter(as.character(.data[[id_col]]) == as.character(.env$phenotype_id))
  if (nrow(selected) != 1L) stop("Expected one phenotype row; found ", nrow(selected), ".")

  sample_columns <- setdiff(names(selected), c(id_columns, id_col, "gene_symbol"))
  if (!is.null(sample_ids)) {
    sample_ids <- as.character(sample_ids)
    if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
      stop("sample_ids must be unique and non-missing.")
    }
    sample_columns <- intersect(sample_ids, sample_columns)
  }
  if (!length(sample_columns)) stop("No matching sample columns were found.")

  result <- selected %>%
    dplyr::select(all_of(sample_columns)) %>%
    pivot_longer(everything(), names_to = "sample_id", values_to = ".value_text",
                 values_transform = list(.value_text = as.character)) %>%
    mutate(.value_text = str_trim(.value_text),
           phenotype = suppressWarnings(as.numeric(.value_text)))

  invalid <- !is.na(result$.value_text) &
    !result$.value_text %in% c("", ".", "NA", "NaN") & is.na(result$phenotype)
  if (any(invalid)) {
    stop("Nonnumeric phenotype values were found. Check the values or supply sample_ids to exclude metadata.")
  }
  n_removed <- sum(!is.finite(result$phenotype))
  result <- result %>% filter(is.finite(phenotype)) %>% dplyr::select(sample_id, phenotype)
  if (!nrow(result)) stop("No finite phenotype values remain.")
  if (n_removed) message("Removed ", n_removed, " missing or non-finite phenotype values.")
  result
}

.sv_tabix_binary <- function(vcf_file, tabix_path) {
  if (!file.exists(vcf_file) || file.access(vcf_file, 4) != 0L) {
    stop("The VCF is not readable: ", vcf_file)
  }
  binary <- if (file.exists(tabix_path)) normalizePath(tabix_path) else unname(Sys.which(tabix_path))
  if (!nzchar(binary) || file.access(binary, 1) != 0L) stop("Cannot execute tabix: ", tabix_path)
  binary
}

.sv_write_region <- function(vcf_file, region, tabix_path, subset_vcf) {
  binary <- .sv_tabix_binary(vcf_file, tabix_path)
  tabix_log <- tempfile(fileext = ".log")
  on.exit(unlink(tabix_log), add = TRUE)
  message("Extracting ", region)
  status <- system2(binary, args = c("-h", shQuote(vcf_file), shQuote(region)),
                    stdout = subset_vcf, stderr = tabix_log)
  if (status != 0L) {
    stop("Tabix failed: ", paste(readLines(tabix_log, warn = FALSE), collapse = "\n"))
  }
}

extract_sv_region <- function(vcf_file, region, tabix_path = "tabix", allow_empty = FALSE,
                              include_samples = TRUE) {
  subset_vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(subset_vcf), add = TRUE)
  .sv_write_region(vcf_file, region, tabix_path, subset_vcf)
  variants <- data.table::fread(
    subset_vcf, skip = "#CHROM", header = TRUE, sep = "\t",
    colClasses = "character", quote = "", na.strings = NULL, check.names = FALSE,
    select = if (!include_samples) c("#CHROM", "POS", "ID", "REF", "ALT", "INFO")
  ) %>%
    as_tibble() %>%
    rename(CHROM = `#CHROM`) %>%
    filter(ALT != ".", !str_detect(ALT, ","))
  if (!nrow(variants) && !allow_empty) stop("No biallelic SVs were returned for ", region, ".")
  variants
}

# Return the original VCF row (including samples) and its 1-based coordinates.
extract_sv_by_id <- function(vcf_file, variant_id, region = NULL, tabix_path = "tabix") {
  if (!is.character(variant_id) || length(variant_id) != 1L || is.na(variant_id) ||
      variant_id %in% c("", ".") || str_detect(variant_id, "[;[:space:]]")) {
    stop("Supply one VCF variant ID.")
  }
  if (!is.null(region)) {
    if (!is.character(region) || length(region) != 1L || is.na(region) ||
        !nzchar(region) || startsWith(region, "-")) stop("Supply one tabix region.")
    subset_vcf <- tempfile(fileext = ".vcf")
    on.exit(unlink(subset_vcf), add = TRUE)
    .sv_write_region(vcf_file, region, tabix_path, subset_vcf)
    vcf_file <- subset_vcf
  }

  con <- gzfile(vcf_file, open = "rt")
  on.exit(close(con), add = TRUE, after = FALSE)
  header <- hit <- character()
  repeat {
    lines <- readLines(con, n = 64L, warn = FALSE)
    if (!length(lines)) break
    chrom_header <- lines[startsWith(lines, "#CHROM\t")]
    if (length(chrom_header)) header <- substring(chrom_header, 2L)
    records <- lines[!startsWith(lines, "#")]
    ids <- str_match(records, "^[^\t]+\t[^\t]+\t([^\t]+)\t")[, 2L]
    matches <- vapply(strsplit(ids, ";", fixed = TRUE), function(x) variant_id %in% x, logical(1))
    hit <- c(hit, records[matches])
    if (length(hit) > 1L) stop("Multiple VCF records match ID ", variant_id, ".")
  }

  if (length(header) != 1L) stop("Missing or invalid #CHROM header in the VCF.")
  if (!length(hit)) {
    where <- if (is.null(region)) "the VCF" else region
    stop("No VCF record matches ID ", variant_id, " in ", where, ".")
  }
  columns <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  fields <- strsplit(hit, "\t", fixed = TRUE)[[1L]]
  if (length(fields) != length(columns)) stop("The matched VCF row and header have different column counts.")
  record <- tibble::as_tibble_row(stats::setNames(fields, columns))
  coordinates <- .sv_cds_spans(record) %>%
    mutate(INFO_END = .sv_info_end(record$INFO), .after = POS)
  list(coordinates = coordinates, record = record)
}

.sv_type_from_vcf <- function(raw_sv) {
  if (!"ALT" %in% names(raw_sv)) stop("SV records must contain an ALT column.")
  alt <- as.character(raw_sv$ALT)
  sv_type <- str_match(alt, "^<([A-Z][A-Z0-9_]*)(?::[^<>,[:space:]]+)?>$")[, 2L]
  source <- rep(NA_character_, length(sv_type))
  source[!is.na(sv_type)] <- "symbolic ALT"

  if ("INFO" %in% names(raw_sv)) {
    annotated_type <- str_match(
      raw_sv$INFO, "(?:^|;)SVTYPE=([A-Z][A-Z0-9_]*)(?::[^;,[:space:]]+)?(?:;|$)"
    )[, 2L]
    annotated <- which(is.na(sv_type) & !is.na(annotated_type))
    sv_type[annotated] <- annotated_type[annotated]
    source[annotated] <- "INFO/SVTYPE"
  }

  # Infer a simple deletion only when the retained sequence matches REF.
  if ("REF" %in% names(raw_sv)) {
    ref <- as.character(raw_sv$REF)
    candidates <- which(is.na(sv_type) & (nchar(ref) - nchar(alt) >= 50L))
    ref_sequence <- toupper(ref[candidates])
    alt_sequence <- toupper(alt[candidates])
    is_deletion <- str_detect(ref_sequence, "^[ACGTN]+$") &
      str_detect(alt_sequence, "^[ACGTN]+$") &
      (startsWith(ref_sequence, alt_sequence) | endsWith(ref_sequence, alt_sequence))
    deletion_rows <- candidates[which(is_deletion)]
    sv_type[deletion_rows] <- "DEL"
    source[deletion_rows] <- "REF/ALT sequences"
  }
  tibble(sv_type = sv_type, sv_type_source = source)
}

# Prefer symbolic ALT, then INFO/SVTYPE, then evidence of a simple deletion.
mask_sv_types <- function(raw_sv, sv_types = NULL) {
  if (is.null(sv_types)) return(raw_sv)
  if (!is.character(sv_types) || !length(sv_types) || anyNA(sv_types) ||
      any(!str_detect(sv_types, "^[A-Z][A-Z0-9_]*$"))) {
    stop('Supply SV type names such as c("DEL", "DUP", "INV").')
  }
  record_types <- .sv_type_from_vcf(raw_sv)$sv_type
  selected <- raw_sv[record_types %in% sv_types, , drop = FALSE]
  if (!nrow(selected)) {
    stop("No extracted SVs match the requested SV types: ", paste(sv_types, collapse = ", "), ".")
  }
  selected
}

# Return a numeric matrix: sample rows, variant columns.
# The default follows the user's merged-callset convention for ./., .|. calls.
# Partial calls, a single '.', and absent GT fields remain NA.
parse_sv_genotypes <- function(raw_sv, sample_ids, missing_as_reference = TRUE) {
  raw_sv <- as_tibble(raw_sv)
  sample_ids <- as.character(sample_ids)
  if (!nrow(raw_sv)) stop("No SV records were supplied.")
  if (!length(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) {
    stop("Supply unique, non-missing sample IDs.")
  }
  required <- c("ID", "ALT", "FORMAT", sample_ids)
  if (!all(required %in% names(raw_sv))) {
    stop("Missing VCF columns: ", paste(setdiff(required, names(raw_sv)), collapse = ", "))
  }
  if (anyNA(raw_sv$ID) || any(raw_sv$ID %in% c("", ".")) || anyDuplicated(raw_sv$ID)) {
    stop("Variant IDs must be unique and non-missing.")
  }
  stopifnot(is.logical(missing_as_reference), length(missing_as_reference) == 1L,
            !is.na(missing_as_reference))

  genotypes <- raw_sv %>%
    mutate(.variant_row = row_number(),
           .gt_index = map_int(str_split(FORMAT, ":"), ~ match("GT", .x))) %>%
    dplyr::select(.variant_row, .gt_index, FORMAT, all_of(sample_ids)) %>%
    pivot_longer(all_of(sample_ids), names_to = "sample_id", values_to = "sample_fields") %>%
    mutate(
      GT = map2_chr(str_split(sample_fields, ":"), .gt_index, ~ {
        if (is.na(.y) || length(.x) < .y) return(NA_character_)
        .x[[.y]]
      }),
      GT = str_remove(GT, "^[/|]"),
      gt_status = case_when(
        is.na(.gt_index) ~ "GT_absent",
        GT %in% c("./.", ".|.") & missing_as_reference ~ "assumed_reference",
        is.na(GT) | GT %in% c("", ".") ~ "missing",
        str_detect(GT, "^[01][/|][01]$") ~ "called",
        str_detect(GT, "^[01.][/|][01.]$") & str_detect(GT, fixed(".")) ~ "missing",
        TRUE ~ "unsupported_GT"
      ),
      dosage = case_when(
        gt_status == "assumed_reference" ~ 0L,
        gt_status == "called" ~ str_count(GT, "1"),
        TRUE ~ NA_integer_
      )
    )
  if (any(genotypes$gt_status == "unsupported_GT")) {
    bad <- genotypes %>% filter(gt_status == "unsupported_GT") %>% distinct(GT) %>% pull(GT)
    stop("Unsupported diploid biallelic GT values: ", paste(head(bad, 5), collapse = ", "))
  }

  Z <- genotypes %>%
    dplyr::select(.variant_row, sample_id, dosage) %>%
    pivot_wider(names_from = sample_id, values_from = dosage) %>%
    arrange(.variant_row) %>%
    dplyr::select(all_of(sample_ids)) %>%
    as.matrix() %>% t()
  colnames(Z) <- raw_sv$ID
  stopifnot(identical(rownames(Z), sample_ids))
  Z
}

# Calculate frequencies only after the final analysis samples have been selected.
filter_sv_genotypes <- function(Z, max_alt_af = 0.01, max_missing = 0.05) {
  stopifnot(is.matrix(Z), is.numeric(Z), nrow(Z) > 0L, ncol(Z) > 0L,
            !is.null(colnames(Z)), !anyDuplicated(colnames(Z)),
            all(Z[!is.na(Z)] %in% 0:2),
            length(max_alt_af) == 1L, is.finite(max_alt_af), max_alt_af > 0, max_alt_af <= 1,
            length(max_missing) == 1L, is.finite(max_missing), max_missing >= 0, max_missing <= 1)
  qc <- tibble(
    variant_id = colnames(Z),
    alt_allele_count = colSums(Z, na.rm = TRUE),
    carriers = colSums(Z > 0, na.rm = TRUE),
    alt_af = colMeans(Z, na.rm = TRUE) / 2,
    missing_rate = colMeans(is.na(Z)),
    variable_genotype = map_lgl(seq_len(ncol(Z)), ~ n_distinct(Z[, .x], na.rm = TRUE) > 1L)
  ) %>%
    mutate(keep = is.finite(alt_af) & alt_af > 0 & alt_af < 1 & variable_genotype &
             alt_af <= max_alt_af & missing_rate <= max_missing)
  summary <- qc %>% summarise(
    total_variants = n(), all_genotypes_missing = sum(missing_rate == 1),
    no_alt_alleles = sum(alt_af == 0, na.rm = TRUE),
    above_af_threshold = sum(alt_af > max_alt_af, na.rm = TRUE),
    above_missing_threshold = sum(missing_rate > max_missing),
    variants_retained = sum(keep)
  )
  list(Z = Z[, qc$keep, drop = FALSE], variant_qc = qc, summary = summary)
}

.sv_model_data <- function(pheno, covariates = NULL, covariate_cols = character()) {
  if (!length(covariate_cols)) {
    if (!is.null(covariates)) stop("Supply covariate_cols when supplying a covariate table.")
    return(list(data = pheno, formula = phenotype ~ 1))
  }
  if (is.null(covariates) || anyDuplicated(covariate_cols) ||
      any(covariate_cols %in% c("sample_id", "phenotype")) ||
      any(make.names(covariate_cols) != covariate_cols)) {
    stop("Supply a covariate table and valid, distinct covariate column names.")
  }
  covariates <- as_tibble(covariates)
  if (!all(c("sample_id", covariate_cols) %in% names(covariates))) {
    stop("Required covariate columns are missing.")
  }
  covariates <- covariates %>% mutate(sample_id = as.character(sample_id))
  if (anyNA(covariates$sample_id) || anyDuplicated(covariates$sample_id)) {
    stop("Covariate sample IDs must be unique and non-missing.")
  }
  data <- pheno %>%
    left_join(dplyr::select(covariates, sample_id, all_of(covariate_cols)), by = "sample_id") %>%
    filter(if_all(all_of(covariate_cols), ~ {
      valid <- !is.na(.x)
      if (is.numeric(.x)) valid <- valid & is.finite(.x)
      valid
    }))
  message("Removed ", nrow(pheno) - nrow(data), " samples with missing or non-finite covariates.")
  list(data = data, formula = stats::reformulate(covariate_cols, response = "phenotype"))
}

.sv_prepare_lof_table <- function(lof_table, lof_missing_as_noncarrier = FALSE) {
  stopifnot(is.logical(lof_missing_as_noncarrier), length(lof_missing_as_noncarrier) == 1L,
            !is.na(lof_missing_as_noncarrier))
  if (is.null(lof_table)) {
    if (lof_missing_as_noncarrier) stop("Supply lof_table to treat unlisted pairs as noncarriers.")
    return(NULL)
  }
  required <- c("sample_id", "gene_id", "has_lof_variant")
  if (!is.data.frame(lof_table) || !all(required %in% names(lof_table))) {
    stop("lof_table must have sample_id, gene_id, and has_lof_variant columns.")
  }
  lof <- lof_table %>%
    as_tibble() %>%
    transmute(
      sample_id = .sv_gene_label(sample_id),
      gene_id = .sv_gene_label(.sv_stable_gene_id(gene_id)),
      has_lof_variant = str_to_lower(.sv_gene_label(has_lof_variant))
    )
  if (anyNA(lof$sample_id) || anyNA(lof$gene_id)) stop("LoF sample and gene IDs must be nonblank.")
  if (any(!is.na(lof$has_lof_variant) & !lof$has_lof_variant %in% c("true", "false", "1", "0", "na"))) {
    stop("Use TRUE, FALSE, 1, 0, or a missing value for has_lof_variant.")
  }
  lof <- lof %>%
    mutate(has_lof_variant = case_when(
      has_lof_variant %in% c("true", "1") ~ 1L,
      has_lof_variant %in% c("false", "0") ~ 0L,
      TRUE ~ NA_integer_
    )) %>%
    distinct()
  if (nrow(distinct(lof, sample_id, gene_id)) != nrow(lof)) {
    stop("Conflicting LoF statuses for the same sample and stable gene ID.")
  }
  lof
}

.sv_with_lof_covariate <- function(model_data, lof_status, lof_missing_as_noncarrier) {
  data <- model_data$data
  if ("has_lof_variant" %in% names(data)) {
    stop("Remove has_lof_variant from covariate_cols when supplying lof_table.")
  }
  matched <- match(data$sample_id, lof_status$sample_id)
  status <- lof_status$has_lof_variant[matched]
  if (lof_missing_as_noncarrier) status[is.na(matched)] <- 0L
  model_data$data <- data %>%
    mutate(has_lof_variant = .env$status) %>%
    filter(!is.na(has_lof_variant))
  use_lof <- n_distinct(model_data$data$has_lof_variant) == 2L
  if (use_lof) {
    model_data$formula <- stats::update(model_data$formula, . ~ . + has_lof_variant)
  }
  model_data$lof_qc <- tibble(
    n_lof_missing = sum(is.na(status)),
    n_lof_carriers = sum(status == 1L, na.rm = TRUE),
    lof_covariate_used = use_lof
  )
  message("Removed ", model_data$lof_qc$n_lof_missing, " samples with unknown LoF status.")
  model_data
}

# Main entry point. No files are read until this function is called.
# covariates: optional table with sample_id and the requested covariate_cols.
# tss_locations: annotation table or TSV path. Maps gene_id to gene_id;
# maps gene_name and Assay to gene_name.
# tss_id / tss_id_col: optional annotation lookup override for assay aliases.
# sv_types: optional vector of SV types to include in one test.
run_sv_set_test <- function(phenotype_df, phenotype_id, vcf_file, chromosome = NULL, tss = NULL,
                            id_col = NULL, flank = 1000000L, tabix_path = "tabix",
                            test = c("SKATO", "burden", "SKAT"),
                            max_alt_af = 0.01, max_missing = 0.05,
                            missing_as_reference = TRUE,
                            covariates = NULL, covariate_cols = character(),
                            tss_locations = NULL, tss_id = NULL, tss_id_col = NULL,
                            ignore_gene_version = FALSE, sv_types = NULL) {
  location <- .sv_resolve_location(
    phenotype_df, phenotype_id, chromosome, tss, id_col,
    tss_locations, tss_id, tss_id_col, ignore_gene_version
  )
  chromosome <- location$chromosome
  tss <- location$tss
  if (!requireNamespace("SKAT", quietly = TRUE)) stop("Install the SKAT package before running a test.")
  test <- match.arg(test)
  region <- make_sv_region(chromosome, tss, flank)
  raw_sv <- extract_sv_region(vcf_file, region, tabix_path)
  result <- .sv_fit_set(
    raw_sv, phenotype_df, phenotype_id, id_col, sv_types, region, test,
    max_alt_af, max_missing, missing_as_reference, covariates, covariate_cols,
    n_svs_extracted = nrow(raw_sv)
  )
  result["tss_location"] <- list(location$annotation)
  result$settings <- c(result$settings, list(
    chromosome = chromosome, tss = tss, flank = flank,
    coordinate_source = if (is.null(tss_locations)) "manual" else "tss_locations",
    tss_id = location$query_id, tss_id_col = location$lookup_col,
    ignore_gene_version = ignore_gene_version
  ))
  result
}

.sv_fit_null_model <- function(model_data) {
  check <- stats::lm(model_data$formula, data = model_data$data, na.action = stats::na.fail)
  if (check$rank < length(stats::coef(check)) || stats::df.residual(check) < 2L ||
      !is.finite(stats::deviance(check)) || stats::deviance(check) <= 0) {
    stop("The null model needs a full-rank design, residual variation, and at least two residual degrees of freedom.")
  }
  SKAT::SKAT_Null_Model(model_data$formula, data = model_data$data, out_type = "C")
}

.sv_run_skat <- function(Z, null_model, test, max_missing) {
  test_args <- list(Z = Z, obj = null_model, weights.beta = c(1, 25),
                    impute.method = "fixed", missing_cutoff = max_missing)
  if (test == "SKATO") {
    test_args$method <- "optimal.adj"
  } else {
    test_args$r.corr <- if (test == "burden") 1 else 0
  }
  do.call(SKAT::SKAT, test_args)
}

.sv_fit_set <- function(raw_sv, phenotype_df, phenotype_id, id_col, sv_types, region, test,
                        max_alt_af, max_missing, missing_as_reference, covariates, covariate_cols,
                        n_svs_extracted, lof_status = NULL, lof_missing_as_noncarrier = FALSE) {
  raw_sv <- mask_sv_types(raw_sv, sv_types)
  mask_label <- if (is.null(sv_types)) "all" else paste(sv_types, collapse = "+")
  fixed_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  pheno <- prep_burden_phenotype(
    phenotype_df, phenotype_id, id_col,
    sample_ids = setdiff(names(raw_sv), fixed_columns)
  )
  model_data <- .sv_model_data(pheno, covariates, covariate_cols)
  if (!is.null(lof_status)) {
    model_data <- .sv_with_lof_covariate(model_data, lof_status, lof_missing_as_noncarrier)
  }
  analysis_df <- model_data$data
  if (nrow(analysis_df) < 3L) stop("Fewer than three usable samples remain.")
  message("Matched analysis samples: ", nrow(analysis_df))

  Z_all <- parse_sv_genotypes(raw_sv, analysis_df$sample_id, missing_as_reference)
  filtered <- filter_sv_genotypes(Z_all, max_alt_af, max_missing)
  print(filtered$summary, width = Inf)
  Z_filtered <- filtered$Z
  if (!ncol(Z_filtered)) stop("No SVs pass the filters. Inspect filter_sv_genotypes() output.")
  if (ncol(Z_filtered) == 1L) warning("One SV remains; this is a single-variant test.")
  stopifnot(identical(rownames(Z_filtered), analysis_df$sample_id))

  null_model <- .sv_fit_null_model(model_data)
  fit <- .sv_run_skat(Z_filtered, null_model, test, max_missing)
  retained <- filtered$variant_qc %>% filter(keep)
  summary <- tibble(
    phenotype_id = as.character(phenotype_id), region = region, test = test,
    sv_mask = mask_label, n_samples = nrow(Z_filtered),
    n_svs_extracted = n_svs_extracted, n_svs_in_mask = ncol(Z_all),
    n_svs_tested = fit$param$n.marker.test,
    singleton_variants = sum(retained$alt_allele_count == 1),
    total_alt_alleles = sum(retained$alt_allele_count),
    median_carriers_per_variant = median(retained$carriers),
    unique_carriers = sum(rowSums(Z_filtered > 0, na.rm = TRUE) > 0),
    p_value = fit$p.value
  )
  summary <- bind_cols(summary, model_data$lof_qc)
  list(
    summary = annotate_sv_results(summary, phenotype_df, id_col),
    phenotype = analysis_df, Z_all = Z_all, Z_filtered = Z_filtered,
    variant_qc = filtered$variant_qc, qc_summary = filtered$summary,
    variants = raw_sv %>% dplyr::select(all_of(fixed_columns)), test_result = fit,
    settings = list(test = test,
                    sv_types = sv_types,
                    max_alt_af = max_alt_af, max_missing = max_missing,
                    missing_as_reference = missing_as_reference,
                    null_formula = paste(deparse(model_data$formula), collapse = " "))
  )
}

# Read the four columns used by the CDS overlap test from a plain or gzipped GTF.
read_cds_gtf <- function(gtf_file) {
  if (!is.character(gtf_file) || length(gtf_file) != 1L || is.na(gtf_file) ||
      !file.exists(gtf_file) || file.access(gtf_file, 4L) != 0L) {
    stop("Supply the path to a readable GTF file.")
  }

  con <- gzfile(gtf_file, open = "rt")
  on.exit(close(con), add = TRUE)
  chunks <- list()

  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    cds <- lines[stringr::str_detect(lines, "^[^#\\t]+\\t[^\\t]*\\tCDS\\t")]
    if (!length(cds)) next

    fields <- stringr::str_split_fixed(cds, "\t", 9L)
    gene_id <- stringr::str_match(
      fields[, 9L], '(?:^|;)\\s*gene_id\\s+"([^"]+)"\\s*(?:;|$)'
    )[, 2L]
    if (anyNA(gene_id)) stop("A CDS row is missing a quoted gene_id attribute.")

    start <- suppressWarnings(as.integer(fields[, 4L]))
    end <- suppressWarnings(as.integer(fields[, 5L]))
    if (anyNA(start) || anyNA(end) ||
        any(!stringr::str_detect(fields[, 4L], "^[0-9]+$") |
            !stringr::str_detect(fields[, 5L], "^[0-9]+$") | start < 1L | end < start)) {
      stop("CDS rows must have positive 1-based integer coordinates with end >= start.")
    }

    chunks[[length(chunks) + 1L]] <- tibble::tibble(
      gene_id = gene_id, chrom = fields[, 1L], start = start, end = end
    )
  }

  if (!length(chunks)) stop("The GTF has no CDS rows.")
  dplyr::bind_rows(chunks)
}

.sv_info_end <- function(info) {
  as.numeric(str_match(info, "(?:^|;)END=([0-9]+)(?:;|$)")[, 2L])
}

.sv_cds_spans <- function(raw) {
  record_type <- .sv_type_from_vcf(raw)
  supported <- record_type$sv_type %in% c("DEL", "DUP", "INV")
  symbolic <- record_type$sv_type_source %in% "symbolic ALT"
  pos <- as.integer(raw$POS)
  info_end <- .sv_info_end(raw$INFO)
  start <- end <- rep(NA_integer_, nrow(raw))

  # Symbolic POS is the padding base and END is inclusive.
  with_end <- which(supported & symbolic & info_end > pos & info_end <= .Machine$integer.max)
  start[with_end] <- pos[with_end] + 1L
  end[with_end] <- as.integer(info_end[with_end])

  # For simple literal deletions, remove the shared prefix or suffix.
  ref <- toupper(raw$REF)
  alt <- toupper(raw$ALT)
  literal_del <- which(!symbolic & record_type$sv_type == "DEL" &
                        nchar(ref) > nchar(alt) & str_detect(ref, "^[ACGTN]+$") &
                        str_detect(alt, "^[ACGTN]+$") &
                        pos + as.double(nchar(ref)) - 1 <= .Machine$integer.max)
  prefix <- startsWith(ref[literal_del], alt[literal_del])
  suffix <- !prefix & endsWith(ref[literal_del], alt[literal_del])
  left_padded <- literal_del[prefix]
  right_padded <- literal_del[suffix]
  start[left_padded] <- pos[left_padded] + nchar(alt[left_padded])
  end[left_padded] <- pos[left_padded] + (nchar(ref[left_padded]) - 1L)
  start[right_padded] <- pos[right_padded]
  end[right_padded] <- pos[right_padded] + (nchar(ref[right_padded]) - nchar(alt[right_padded]) - 1L)

  tibble(
    CHROM = raw$CHROM, POS = pos, variant_id = raw$ID, ALT = raw$ALT,
    sv_type = record_type$sv_type, affected_start = start, affected_end = end,
    span_status = case_when(!supported ~ "unsupported_type",
                            is.na(start) ~ "unresolved_span", TRUE ~ "ok")
  )
}

.sv_load_cds <- function(gtf) {
  cds <- if (is.character(gtf) && length(gtf) == 1L) read_cds_gtf(gtf) else as_tibble(gtf)
  if (!all(c("gene_id", "chrom", "start", "end") %in% names(cds))) {
    stop("Use a GTF path or a CDS table with gene_id, chrom, start, and end columns.")
  }
  cds
}

.sv_cds_gene_key <- function(gene_id, ignore_gene_version) {
  if (ignore_gene_version) str_remove(.sv_stable_gene_id(gene_id), "\\.[0-9]+$") else as.character(gene_id)
}

# Query the range covered by a gene's CDS, then require an actual CDS overlap.
extract_cds_svs <- function(vcf_file, gtf, gene_id, tabix_path = "tabix",
                            ignore_gene_version = FALSE) {
  if (!requireNamespace("IRanges", quietly = TRUE)) stop("Install the Bioconductor package IRanges.")
  cds <- .sv_load_cds(gtf)
  stopifnot(length(gene_id) == 1L, !is.na(gene_id), nzchar(gene_id),
            is.logical(ignore_gene_version), length(ignore_gene_version) == 1L, !is.na(ignore_gene_version))
  query_key <- .sv_cds_gene_key(gene_id, ignore_gene_version)
  selected <- cds %>%
    dplyr::filter(.sv_cds_gene_key(gene_id, .env$ignore_gene_version) == .env$query_key)
  if (!nrow(selected)) stop("No CDS rows match gene ID ", gene_id, ".")
  chrom <- unique(as.character(selected$chrom))
  if (n_distinct(selected$gene_id) != 1L || length(chrom) != 1L ||
      is.na(chrom) || !nzchar(chrom) || str_detect(chrom, "[:\\s]")) {
    stop("The CDS lookup must resolve to one gene ID and chromosome.")
  }
  start <- selected$start
  end <- selected$end
  if (!is.numeric(start) || !is.numeric(end) ||
      any(!is.finite(start) | !is.finite(end) | start != floor(start) | end != floor(end) |
          start < 1 | end < start | end > .Machine$integer.max)) {
    stop("CDS coordinates must be positive integers with end >= start.")
  }
  union_cds <- IRanges::reduce(IRanges::IRanges(start, end))
  bounds <- sprintf("%s:%d-%d", chrom, min(start), max(end))
  raw <- extract_sv_region(vcf_file, bounds, tabix_path, allow_empty = TRUE)
  overlap <- .sv_cds_spans(raw)
  known <- which(overlap$span_status == "ok" & overlap$CHROM == chrom)
  overlap$cds_overlap_bp <- rep(NA_real_, nrow(overlap))
  overlap$cds_overlap_bp[known] <- 0
  if (length(known)) {
    spans <- IRanges::IRanges(overlap$affected_start[known], overlap$affected_end[known])
    hits <- IRanges::findOverlaps(spans, union_cds)
    q <- S4Vectors::queryHits(hits)
    s <- S4Vectors::subjectHits(hits)
    if (length(q)) {
      bp <- pmin(IRanges::end(spans)[q], IRanges::end(union_cds)[s]) -
        pmax(IRanges::start(spans)[q], IRanges::start(union_cds)[s]) + 1
      totals <- rowsum(bp, q)
      overlap$cds_overlap_bp[known[as.integer(rownames(totals))]] <- as.numeric(totals)
    }
  }
  list(
    region = bounds, gene_id = selected$gene_id[[1L]],
    cds_intervals = tibble(chrom = chrom, start = IRanges::start(union_cds), end = IRanges::end(union_cds)),
    overlap = overlap, variants = raw[which(overlap$cds_overlap_bp > 0), , drop = FALSE]
  )
}

run_sv_cds_test <- function(phenotype_df, phenotype_id, vcf_file, gtf, cds_gene_id = NULL,
                            id_col = "gene_id", tabix_path = "tabix", ignore_gene_version = FALSE,
                            test = c("SKATO", "burden", "SKAT"), sv_types = c("DEL", "DUP"),
                            max_alt_af = 0.01, max_missing = 0.05, missing_as_reference = TRUE,
                            covariates = NULL, covariate_cols = character(),
                            lof_table = NULL, lof_missing_as_noncarrier = FALSE) {
  if (!requireNamespace("SKAT", quietly = TRUE)) stop("Install the SKAT package before running a test.")
  test <- match.arg(test)
  if (!is.null(sv_types) && (!length(sv_types) || any(!sv_types %in% c("DEL", "DUP", "INV")))) {
    stop("CDS interval selection supports DEL, DUP, and INV.")
  }
  if (is.null(cds_gene_id)) {
    cds_gene_id <- .sv_harmonized_gene_id(phenotype_df, phenotype_id, id_col)
    if (is.null(cds_gene_id)) cds_gene_id <- phenotype_id
  }
  lof_status <- .sv_prepare_lof_table(lof_table, lof_missing_as_noncarrier)
  if (!is.null(lof_status)) {
    lof_gene_id <- .sv_stable_gene_id(cds_gene_id)
    lof_status <- lof_status %>% filter(gene_id == .env$lof_gene_id)
  }
  cds <- extract_cds_svs(vcf_file, gtf, cds_gene_id, tabix_path, ignore_gene_version)
  if (!nrow(cds$variants)) {
    stop("No SVs with resolvable spans overlap the CDS. Inspect extract_cds_svs() output.")
  }
  result <- .sv_fit_set(
    cds$variants, phenotype_df, phenotype_id, id_col, sv_types, cds$region, test,
    max_alt_af, max_missing, missing_as_reference, covariates, covariate_cols,
    n_svs_extracted = nrow(cds$overlap),
    lof_status = lof_status, lof_missing_as_noncarrier = lof_missing_as_noncarrier
  )
  active_types <- if (is.null(sv_types)) c("DEL", "DUP", "INV") else sv_types
  result$summary$n_svs_overlapping_cds <- nrow(cds$variants)
  result$summary$n_svs_unresolved_span <- sum(cds$overlap$span_status == "unresolved_span" &
                                              cds$overlap$sv_type %in% active_types)
  result$cds_intervals <- cds$cds_intervals
  result$cds_overlap <- cds$overlap
  result$settings$coordinate_source <- "gtf_cds"
  result$settings$cds_gene_id <- cds$gene_id
  result$settings$ignore_gene_version <- ignore_gene_version
  if (!is.null(lof_status)) result$settings$lof_missing_as_noncarrier <- lof_missing_as_noncarrier
  result
}

# Run one CDS test per requested gene. Use harmonize_sv_inputs() to add gene symbols.
run_sv_cds_tests <- function(phenotype_df, gene_ids, vcf_file, gtf, ...,
                             id_col = "gene_id", tabix_path = "tabix", ignore_gene_version = FALSE,
                             test = c("SKATO", "burden", "SKAT"), follow_up = FALSE, q_cutoff = 0.05,
                             keep_carrier_scores = FALSE, original_vcf_file = NULL,
                             lof_table = NULL, lof_missing_as_noncarrier = FALSE) {
  phenotype_df <- as_tibble(phenotype_df)
  if (!"gene_id" %in% names(phenotype_df)) stop("The phenotype table needs a gene_id column.")
  if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(phenotype_df)) {
    stop("id_col must name one existing phenotype column.")
  }
  if (!is.character(gene_ids) || !length(gene_ids) || anyNA(.sv_gene_label(gene_ids))) {
    stop("Supply a nonempty character vector of gene IDs.")
  }
  stopifnot(is.logical(ignore_gene_version), length(ignore_gene_version) == 1L, !is.na(ignore_gene_version))
  stopifnot(is.logical(follow_up), length(follow_up) == 1L, !is.na(follow_up))
  stopifnot(is.logical(keep_carrier_scores), length(keep_carrier_scores) == 1L, !is.na(keep_carrier_scores))
  if (!is.numeric(q_cutoff) || length(q_cutoff) != 1L || !is.finite(q_cutoff) ||
      q_cutoff < 0 || q_cutoff > 1) stop("q_cutoff must be one number from 0 to 1.")
  ids <- unique(.sv_stable_gene_id(gene_ids))
  shared_options <- list(...)
  allowed_options <- setdiff(names(formals(run_sv_cds_test)), c(
    "phenotype_df", "phenotype_id", "vcf_file", "gtf", "cds_gene_id", "id_col",
    "tabix_path", "ignore_gene_version", "test", "lof_table", "lof_missing_as_noncarrier"
  ))
  if (length(shared_options) && (is.null(names(shared_options)) ||
      any(!names(shared_options) %in% allowed_options))) {
    stop("Name shared CDS test options, such as sv_types or max_alt_af.")
  }
  if (!requireNamespace("SKAT", quietly = TRUE)) stop("Install the SKAT package before running tests.")
  if (!requireNamespace("qvalue", quietly = TRUE)) stop("Install the Bioconductor package qvalue before running CDS tests.")
  tabix_path <- .sv_tabix_binary(vcf_file, tabix_path)
  test <- match.arg(test)
  cds <- .sv_load_cds(gtf)
  lof <- .sv_prepare_lof_table(lof_table, lof_missing_as_noncarrier)
  lof_rows <- if (!is.null(lof)) split(seq_len(nrow(lof)), lof$gene_id)

  phenotype_rows <- split(seq_len(nrow(phenotype_df)), factor(.sv_stable_gene_id(phenotype_df$gene_id), levels = ids))
  cds_rows <- split(seq_len(nrow(cds)), .sv_cds_gene_key(cds$gene_id, ignore_gene_version))
  symbols <- purrr::map_chr(ids, function(id) {
    labels <- .sv_gene_label(phenotype_df[["gene_symbol"]][phenotype_rows[[id]]])
    if (all(is.na(labels))) labels <- .sv_gene_label(cds[["gene_symbol"]][cds_rows[[id]]])
    labels <- unique(labels[!is.na(labels)])
    if (length(labels) == 1L) labels[[1L]] else NA_character_
  })
  save_result <- NULL
  if (keep_carrier_scores) carrier_scores <- set_names(vector("list", length(ids)), ids)
  if (follow_up) {
    scratch <- tempfile("sv-cds-followup-")
    if (!dir.create(scratch, mode = "0700")) stop("Cannot create temporary storage for the CDS follow-up.")
    on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
    fit_paths <- file.path(scratch, paste0(seq_along(ids), ".rds"))
  }
  if (follow_up || keep_carrier_scores) {
    save_result <- function(result, index) {
      if (follow_up) {
        inputs <- result[c("summary", "phenotype", "Z_filtered", "variants", "settings")]
        inputs$test_result <- result$test_result["test.snp.mac"]
        saveRDS(inputs, fit_paths[[index]])
      }
      if (keep_carrier_scores) {
        carrier_scores[[index]] <<- tryCatch(
          extract_sv_carrier_zscores(result) %>%
            mutate(gene_id = .env$ids[[index]], gene_symbol = .env$symbols[[index]], .after = phenotype_id),
          error = conditionMessage
        )
      }
    }
  }
  results <- .sv_collect_batch(ids, test, function(id) {
    rows <- phenotype_rows[[id]]
    if (!length(rows)) stop("No phenotype row for gene ID ", id, ".")
    if (length(rows) != 1L) stop("Found multiple phenotype rows for gene ID ", id, ".")
    phenotype_id <- as.character(phenotype_df[[id_col]][rows])
    if (is.na(.sv_gene_label(phenotype_id))) stop("Missing phenotype identifier for gene ID ", id, ".")
    cds_gene_id <- as.character(phenotype_df$gene_id[rows])
    cds_key <- .sv_cds_gene_key(cds_gene_id, ignore_gene_version)
    run_sv_cds_test(
      phenotype_df = phenotype_df[rows, , drop = FALSE], phenotype_id = phenotype_id,
      vcf_file = vcf_file, gtf = cds[cds_rows[[cds_key]], , drop = FALSE],
      cds_gene_id = cds_gene_id, id_col = id_col, tabix_path = tabix_path,
      ignore_gene_version = ignore_gene_version, test = test,
      lof_table = if (!is.null(lof)) lof[lof_rows[[id]], , drop = FALSE],
      lof_missing_as_noncarrier = lof_missing_as_noncarrier, ...
    )
  }, progress = "SV CDS tests", save_result = save_result)

  tested <- which(results$status == "ok" & results[["n_svs_tested"]] > 0L)
  results$q_value <- NA_real_
  if (length(tested)) results$q_value[tested] <- .sv_storey_qvalues(results$p_value[tested])

  results <- results %>%
    mutate(gene_id = .env$ids, gene_symbol = .env$symbols) %>%
    relocate(gene_id, gene_symbol)
  if (!follow_up && !keep_carrier_scores) return(results)
  batch <- list(summary = results)
  if (keep_carrier_scores) batch$carrier_scores <- carrier_scores
  if (!follow_up) return(batch)

  selected <- which(results$q_value <= q_cutoff)
  followups <- purrr::map(selected, function(index) {
    tryCatch(
      run_sv_set_followup(readRDS(fit_paths[[index]])),
      error = function(e) {
        tibble(phenotype_id = results$phenotype_id[[index]], variant_id = NA_character_) %>%
          .sv_target_output() %>% mutate(error_message = conditionMessage(e))
      }
    ) %>% mutate(gene_id = .env$ids[[index]], gene_symbol = .env$symbols[[index]])
  }, .progress = "CDS single-variant follow-up")
  empty <- tibble(gene_id = character(), gene_symbol = character(), phenotype_id = character(),
                  variant_id = character()) %>% .sv_target_output()
  batch$single_variant <- bind_rows(empty, followups) %>% dplyr::select(-region)
  if (!is.null(original_vcf_file)) {
    batch$single_variant <- restore_sv_followup_ids(batch$single_variant, original_vcf_file, tabix_path)
  }
  batch
}

.sv_storey_qvalues <- function(p_values) {
  # Keep the default grid points supported by the observed p-values.
  lambda <- seq(0.05, 0.95, 0.05)
  lambda <- lambda[lambda <= max(p_values)]
  tryCatch({
    if (length(lambda) < 4L) stop("Too few supported points on the default lambda grid.")
    qvalue::qvalue(p_values, lambda = lambda, lfdr.out = FALSE)$qvalues
  }, error = function(e) {
    warning("Storey q-values are unavailable. ", conditionMessage(e), call. = FALSE)
    rep(NA_real_, length(p_values))
  })
}

# Run one test per identifier. Pass shared test options as named arguments.
# Keep summary rows only, and adjust p-values across successful tests.
run_sv_set_tests <- function(phenotype_df, vcf_file, tss_locations, ...,
                             id_col = "gene_id", tabix_path = "tabix",
                             test = c("SKATO", "burden", "SKAT")) {
  phenotype_df <- as_tibble(phenotype_df)
  if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(phenotype_df)) {
    stop("id_col must name one existing phenotype column.")
  }
  ids <- as.character(phenotype_df[[id_col]])
  if (!length(ids) || anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop("The phenotype table must contain non-missing, non-blank identifiers.")
  }
  groups <- split(seq_along(ids), factor(ids, levels = unique(ids)))
  shared_args <- list(...)
  if (length(shared_args) &&
      (is.null(names(shared_args)) || any(!nzchar(names(shared_args))))) {
    stop("Name shared test options, such as flank or sv_types.")
  }
  allowed_args <- setdiff(names(formals(run_sv_set_test)), c(
    "phenotype_df", "phenotype_id", "vcf_file", "chromosome", "tss", "id_col",
    "tabix_path", "test", "tss_locations", "tss_id"
  ))
  unsupported_args <- setdiff(names(shared_args), allowed_args)
  if (length(unsupported_args)) {
    stop("Unsupported batch options: ", paste(unsupported_args, collapse = ", "), ".")
  }
  if (is.character(tss_locations) && length(tss_locations) == 1L) {
    tss_locations <- data.table::fread(file = tss_locations, colClasses = "character", check.names = FALSE)
  }
  if (!is.data.frame(tss_locations)) stop("Supply an annotation table or a file path as tss_locations.")
  if (!requireNamespace("SKAT", quietly = TRUE)) stop("Install the SKAT package before running tests.")
  tabix_path <- .sv_tabix_binary(vcf_file, tabix_path)
  test <- match.arg(test)

  results <- .sv_collect_batch(names(groups), test, function(id) {
    rows <- groups[[id]]
    if (length(rows) != 1L) stop("Found multiple phenotype rows for ", id, ".")
    run_sv_set_test(
      phenotype_df = phenotype_df[rows, , drop = FALSE], phenotype_id = id,
      vcf_file = vcf_file, tss_locations = tss_locations, id_col = id_col,
      tabix_path = tabix_path, test = test, ...
    )
  }, progress = "SV set tests")
  tested <- results$status == "ok"
  results$q_value <- NA_real_
  results$q_value[tested] <- stats::p.adjust(results$p_value[tested], method = "BH")
  annotate_sv_results(results, phenotype_df, id_col)
}

.sv_collect_batch <- function(ids, test, run_one, progress, save_result = NULL) {
  entries <- purrr::map(seq_along(ids), function(index) {
    id <- ids[[index]]
    warnings <- character()
    summary <- tryCatch(
      withCallingHandlers({
        invisible(utils::capture.output({
          result <- suppressMessages(run_one(id))
        }))
        p_value <- result$summary$p_value
        if (length(p_value) != 1L || !is.finite(p_value) || p_value < 0 || p_value > 1) {
          stop("The test did not return a valid p-value.")
        }
        result$summary %>% mutate(status = "ok", error_message = NA_character_)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) tibble(
        phenotype_id = id, test = test, p_value = NA_real_,
        status = "error", error_message = conditionMessage(e)
      )
    )
    summary$warnings <- if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_
    if (!is.null(save_result) && summary$status == "ok") save_result(result, index)
    summary
  }, .progress = progress)
  bind_rows(entries)
}

.sv_target_output <- function(rows) {
  rows %>% mutate(
    region = NA_character_, CHROM = NA_character_, POS = NA_integer_, ALT = NA_character_,
    n_samples = NA_integer_, carriers = NA_integer_, alt_af = NA_real_, missing_rate = NA_real_,
    beta = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p_value = NA_real_,
    status = "error", error_message = NA_character_
  )
}

.sv_tested_variant_ids <- function(set_result) {
  summary <- set_result[["summary"]]
  n_tested <- summary[["n_svs_tested"]]
  if (!is.data.frame(summary) || nrow(summary) != 1L || length(n_tested) != 1L ||
      !is.finite(n_tested) || n_tested < 0L) {
    stop("Supply the full result from one run_sv_set_test() or run_sv_cds_test() call.")
  }

  # SKAT can exclude columns that passed our filters. Its named output gives the actual membership.
  skat_mac <- set_result$test_result[["test.snp.mac"]]
  ids <- as.character(names(skat_mac))
  Z <- set_result$Z_filtered

  # SKAT can drop the column name when only one SV remains.
  if (!length(ids) && n_tested == 1L && length(skat_mac) == 1L) {
    candidates <- colnames(Z)
    if (length(candidates) > 1L) {
      candidates <- candidates[colMeans(is.na(Z)) < set_result$settings$max_missing]
    }
    ids <- candidates
  }
  if (length(ids) != n_tested || anyNA(ids) || anyDuplicated(ids) ||
      any(!ids %in% colnames(Z)) || any(!ids %in% set_result$variants$ID)) {
    stop("The SKAT variant IDs do not match the saved genotypes and variant records.")
  }
  ids
}

# Pool the saved expression z-scores across the tested SVs of each class.
# Each sample contributes once per class, using observed alternate-allele calls.
# Keep the VCF IDs on each row as a comma-delimited string.
# For a CDS batch or its RDS file, pass one gene_id.
extract_sv_carrier_zscores <- function(set_result, gene_id = NULL) {
  if (is.character(set_result) && length(set_result) == 1L) {
    set_result <- tryCatch(readRDS(set_result), error = function(e) {
      stop("Cannot read the result as an RDS file. ", conditionMessage(e), call. = FALSE)
    })
  }
  if (!is.list(set_result) || is.null(set_result[["Z_filtered"]])) {
    scores <- if (is.list(set_result)) set_result[["carrier_scores"]]
    if (is.null(scores)) {
      stop("This result has no carrier scores. Use run_sv_cds_tests() with keep_carrier_scores = TRUE.")
    }
    if (!is.character(gene_id) || length(gene_id) != 1L || is.na(.sv_gene_label(gene_id))) {
      stop("Supply one Ensembl gene_id from the CDS batch.")
    }
    id <- .sv_stable_gene_id(gene_id)
    if (!id %in% names(scores)) stop("Gene ", id, " is absent from the CDS batch.")
    carriers <- scores[[id]]
    if (is.null(carriers)) stop("Gene ", id, " has no saved carrier scores. Check its batch summary.")
    if (is.character(carriers)) stop("Cannot extract carrier scores for ", id, ". ", carriers)
    return(carriers)
  }
  if (!is.null(gene_id)) stop("Use gene_id with a CDS batch result.")

  ids <- .sv_tested_variant_ids(set_result)
  Z <- set_result$Z_filtered
  data <- set_result$phenotype
  if (!identical(rownames(Z), data$sample_id)) stop("The saved sample and genotype orders do not match.")
  if (!is.numeric(data$phenotype) || any(!is.finite(data$phenotype))) {
    stop("The saved phenotype must contain finite numeric expression z-scores.")
  }

  variants <- set_result$variants[match(ids, set_result$variants$ID), , drop = FALSE]
  sv_types <- .sv_type_from_vcf(variants)$sv_type
  gene <- set_result$summary %>% dplyr::select(any_of(c("phenotype_id", "gene_id", "gene_symbol")))
  carriers <- map(c("DEL", "DUP"), function(type) {
    type_ids <- ids[which(sv_types == type)]
    has_alt <- Z[, type_ids, drop = FALSE] > 0L
    counts <- rowSums(has_alt, na.rm = TRUE)
    rows <- which(unname(counts) > 0L)
    tibble(
      sample_id = data$sample_id[rows], sv_type = type, z_score = data$phenotype[rows],
      n_carried_svs = as.integer(counts[rows]),
      variant_ids = map_chr(rows, function(row) {
        str_c(type_ids[which(has_alt[row, ])], collapse = ",")
      })
    )
  }) %>% list_rbind()
  bind_cols(gene[rep(1L, nrow(carriers)), , drop = FALSE], carriers)
}

# Fit each SV that SKAT used, with the original samples and covariates.
run_sv_set_followup <- function(set_result, original_vcf_file = NULL, tabix_path = "tabix") {
  ids <- .sv_tested_variant_ids(set_result)
  summary <- set_result$summary
  n_tested <- summary[["n_svs_tested"]]
  Z <- set_result$Z_filtered
  variants <- set_result$variants[match(ids, set_result$variants$ID), , drop = FALSE]
  data <- set_result$phenotype
  if (!identical(rownames(Z), data$sample_id)) stop("The saved sample and genotype orders do not match.")
  out <- summary[rep(1L, n_tested), , drop = FALSE] %>%
    dplyr::select(any_of(c("phenotype_id", "gene_id", "gene_symbol"))) %>%
    mutate(variant_id = as.character(.env$variants$ID)) %>%
    .sv_target_output() %>% dplyr::select(-region) %>%
    mutate(CHROM = variants$CHROM, POS = as.integer(variants$POS), ALT = variants$ALT)
  if ("REF" %in% names(variants)) out$REF <- as.character(variants$REF)
  if ("INFO" %in% names(variants)) out$INFO_END <- .sv_info_end(variants$INFO)
  if (!is.null(original_vcf_file)) out <- restore_sv_followup_ids(out, original_vcf_file, tabix_path)
  if (!length(ids)) return(out)

  formula <- stats::as.formula(set_result$settings$null_formula)
  if (".sv_dosage" %in% all.vars(formula)) stop("Choose a different name for the .sv_dosage covariate.")
  model_data <- list(data = data, formula = stats::update(formula, . ~ . + .sv_dosage))
  dosage <- Z[, ids, drop = FALSE]
  out <- out %>% mutate(
    n_samples = nrow(dosage), carriers = as.integer(colSums(dosage > 0, na.rm = TRUE)),
    alt_af = colMeans(dosage, na.rm = TRUE) / 2, missing_rate = colMeans(is.na(dosage))
  )
  purrr::map(seq_along(ids), function(i) {
    row <- out[i, ]
    tryCatch(.sv_fit_target_dosage(row, dosage[, i], model_data), error = function(e) {
      row$error_message <- conditionMessage(e)
      row
    })
  }) %>% purrr::list_rbind()
}

# Recover the ID from the original VCF for each saved follow-up row.
restore_sv_followup_ids <- function(followup, vcf_file, tabix_path = "tabix") {
  out <- as_tibble(followup) %>% ungroup()
  required <- c("CHROM", "POS", "ALT")
  if (!all(required %in% names(out))) stop("The follow-up table needs CHROM, POS, and ALT.")
  old_ids <- as.character(out[["variant_id"]])
  if (!length(old_ids)) old_ids <- rep(NA_character_, nrow(out))
  if ("previous_variant_id" %in% names(out)) {
    old_ids <- dplyr::coalesce(old_ids, as.character(out$previous_variant_id))
  }
  out <- out %>% dplyr::select(-any_of(c("region", "previous_variant_id", "skat_variant_id")))
  out$variant_id <- rep(NA_character_, nrow(out))
  out$id_match_status <- rep("missing_coordinates", nrow(out))
  loci <- out %>%
    transmute(row = row_number(), chrom = as.character(CHROM), pos = as.numeric(as.character(POS)), ALT) %>%
    filter(!is.na(chrom), nzchar(chrom), is.finite(pos), pos >= 1, pos == floor(pos), !is.na(ALT), ALT != ".") %>%
    group_by(chrom, pos) %>% group_split()

  for (locus in loci) {
    chrom <- locus$chrom[[1L]]
    pos <- locus$pos[[1L]]
    raw <- extract_sv_region(vcf_file, make_sv_region(chrom, pos, flank = 0L), tabix_path,
                             allow_empty = TRUE, include_samples = FALSE) %>%
      filter(CHROM == .env$chrom, as.numeric(POS) == .env$pos)
    for (i in locus$row) {
      candidates <- raw[which(raw$ALT == out$ALT[[i]]), , drop = FALSE]
      if ("REF" %in% names(out) && !is.na(out$REF[[i]]) && out$REF[[i]] != ".") {
        candidates <- candidates[which(candidates$REF == out$REF[[i]]), , drop = FALSE]
      }
      if ("INFO_END" %in% names(out) && !is.na(out$INFO_END[[i]])) {
        candidates <- candidates[which(.sv_info_end(candidates$INFO) == out$INFO_END[[i]]), , drop = FALSE]
      }
      vcf_ids <- as.character(candidates$ID)
      if (length(vcf_ids) > 1L) {
        known_id <- old_ids[i]
        if (!is.na(known_id) && !known_id %in% c("", ".")) {
          same_id <- vcf_ids %in% known_id | map_lgl(strsplit(vcf_ids, ";", fixed = TRUE), ~ known_id %in% .x)
          if (any(same_id)) vcf_ids <- vcf_ids[same_id]
        }
      }
      if (length(vcf_ids) == 1L && !is.na(vcf_ids) && !vcf_ids %in% c("", ".")) {
        out$variant_id[i] <- vcf_ids
        out$id_match_status[i] <- "matched"
      } else {
        out$id_match_status[i] <- case_when(
          length(vcf_ids) == 0L ~ "not_found", length(vcf_ids) > 1L ~ "ambiguous", TRUE ~ "missing_vcf_id"
        )
      }
    }
  }
  out
}

.sv_target_model <- function(phenotype_df, phenotype_id, id_col, raw, covariates, covariate_cols) {
  fixed_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  if (".sv_dosage" %in% covariate_cols) stop("Choose a different name for the .sv_dosage covariate.")
  pheno <- prep_burden_phenotype(
    phenotype_df, phenotype_id, id_col, sample_ids = setdiff(names(raw), fixed_columns)
  )
  model_data <- .sv_model_data(pheno, covariates, covariate_cols)
  if (nrow(model_data$data) < 3L) stop("Fewer than three usable samples remain.")
  model_data$formula <- stats::update(model_data$formula, . ~ . + .sv_dosage)
  model_data
}

.sv_fit_target_record <- function(row, record, model_data, max_alt_af, max_missing, missing_as_reference) {
  tryCatch({
    if (nrow(record) != 1L) stop("The VCF contains multiple records with this ID in the window.")
    row$CHROM <- record$CHROM[[1L]]
    row$POS <- as.integer(record$POS[[1L]])
    row$ALT <- record$ALT[[1L]]
    row$n_samples <- nrow(model_data$data)

    # The parser needs a column label even when this one VCF record has no ID.
    if (is.na(record$ID[[1L]]) || record$ID[[1L]] %in% c("", ".")) record$ID <- "target_sv"
    Z <- parse_sv_genotypes(record, model_data$data$sample_id, missing_as_reference)
    qc <- filter_sv_genotypes(Z, max_alt_af, max_missing)$variant_qc
    dosage <- Z[, 1L]
    row$carriers <- sum(dosage > 0, na.rm = TRUE)
    row$alt_af <- unname(qc$alt_af)
    row$missing_rate <- unname(qc$missing_rate)
    if (!qc$keep) {
      row$status <- "filtered"
      row$error_message <- "The SV failed the frequency, genotype variation, or missingness filter."
      return(row)
    }

    .sv_fit_target_dosage(row, dosage, model_data)
  }, error = function(e) {
    row$error_message <- conditionMessage(e)
    row
  })
}

.sv_fit_target_dosage <- function(row, dosage, model_data) {
  dat <- model_data$data %>%
    mutate(.sv_dosage = coalesce(.env$dosage, mean(.env$dosage, na.rm = TRUE)))
  fit <- stats::lm(model_data$formula, data = dat, na.action = stats::na.fail)
  if (fit$rank < length(stats::coef(fit)) || stats::df.residual(fit) < 1L ||
      stats::deviance(fit) <= 0) stop("The model cannot estimate the SV coefficient.")
  estimate <- broom::tidy(fit, conf.int = TRUE) %>% dplyr::filter(term == ".sv_dosage")
  if (nrow(estimate) != 1L || !is.finite(estimate$p.value)) {
    stop("The model did not return a finite SV p-value.")
  }
  row$beta <- estimate$estimate
  row$conf.low <- estimate$conf.low
  row$conf.high <- estimate$conf.high
  row$p_value <- estimate$p.value
  row$status <- "ok"
  row
}

# Test requested gene and SV pairs using exact VCF IDs in each gene window.
# max_alt_af = 1 lets the requested IDs define the frequency range.
run_targeted_sv_tests <- function(targets, phenotype_df, vcf_file, tss_locations,
                                  id_col = "gene_id", flank = 100000L, tabix_path = "tabix",
                                  max_alt_af = 1, max_missing = 0.05,
                                  missing_as_reference = TRUE, covariates = NULL,
                                  covariate_cols = character(), tss_id_col = NULL,
                                  ignore_gene_version = FALSE) {
  targets <- as_tibble(targets)
  if (!all(c("phenotype_id", "variant_id") %in% names(targets))) {
    stop("Targets must have phenotype_id and variant_id columns.")
  }
  targets <- targets %>%
    transmute(phenotype_id = as.character(phenotype_id), variant_id = as.character(variant_id)) %>%
    distinct()
  if (!nrow(targets) || anyNA(targets) ||
      any(!nzchar(trimws(targets$phenotype_id))) || any(!nzchar(trimws(targets$variant_id)))) {
    stop("Supply at least one pair of nonblank phenotype and variant IDs.")
  }
  if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(phenotype_df)) {
    stop("id_col must name one existing phenotype column.")
  }
  if (is.character(tss_locations) && length(tss_locations) == 1L) {
    tss_locations <- data.table::fread(file = tss_locations, colClasses = "character", check.names = FALSE)
  }
  tabix_path <- .sv_tabix_binary(vcf_file, tabix_path)
  by_gene <- split(seq_len(nrow(targets)), factor(targets$phenotype_id, levels = unique(targets$phenotype_id)))

  purrr::map(by_gene, function(rows) {
    out <- targets[rows, ] %>% mutate(.input_order = .env$rows) %>% .sv_target_output()
    tryCatch({
      id <- out$phenotype_id[[1L]]
      location <- .sv_resolve_location(
        phenotype_df, id, id_col = id_col, tss_locations = tss_locations,
        tss_id_col = tss_id_col, ignore_gene_version = ignore_gene_version
      )
      out$region <- rep(make_sv_region(location$chromosome, location$tss, flank), nrow(out))
      raw <- extract_sv_region(vcf_file, out$region[[1L]], tabix_path, allow_empty = TRUE)
      in_window <- out$variant_id %in% raw$ID
      out$status[!in_window] <- "not_in_window"
      if (!any(in_window)) return(out)

      model_data <- .sv_target_model(phenotype_df, id, id_col, raw, covariates, covariate_cols)

      purrr::map(seq_len(nrow(out)), function(i) {
        row <- out[i, ]
        if (!in_window[[i]]) return(row)
        record <- raw[which(raw$ID == row$variant_id[[1L]]), , drop = FALSE]
        .sv_fit_target_record(row, record, model_data, max_alt_af, max_missing, missing_as_reference)
      }) %>% purrr::list_rbind()
    }, error = function(e) {
      out$error_message[out$status == "error"] <- conditionMessage(e)
      out
    })
  }, .progress = "Targeted SV tests") %>%
    purrr::list_rbind() %>% arrange(.input_order) %>% dplyr::select(-.input_order) %>%
    annotate_sv_results(phenotype_df, id_col)
}

# Test the phenotype and exact VCF position specified on each sheet row.
run_sv_sheet_tests <- function(sv_sheet, phenotype_df, vcf_file,
                               gene_col = "gene_id", chrom_col = "chrom", pos_col = "pos",
                               variant_id_col = NULL, id_col = "gene_id", tabix_path = "tabix",
                               max_alt_af = 1, max_missing = 0.05, missing_as_reference = TRUE,
                               covariates = NULL, covariate_cols = character()) {
  sv_sheet <- as_tibble(sv_sheet)
  required <- c(gene_col, chrom_col, pos_col, variant_id_col)
  absent <- setdiff(required, names(sv_sheet))
  if (length(absent)) stop("Missing sheet columns: ", paste(absent, collapse = ", "), ".")
  if (!nrow(sv_sheet)) stop("The SV sheet has no rows.")
  if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(phenotype_df)) {
    stop("id_col must name one existing phenotype column.")
  }
  tabix_path <- .sv_tabix_binary(vcf_file, tabix_path)

  requested_ids <- if (is.null(variant_id_col)) NA_character_ else trimws(as.character(sv_sheet[[variant_id_col]]))
  requested_ids[requested_ids %in% c("", ".")] <- NA_character_
  targets <- sv_sheet %>% transmute(
    sheet_row = row_number(), phenotype_id = as.character(.data[[gene_col]]),
    query_chrom = trimws(as.character(.data[[chrom_col]])),
    query_pos = suppressWarnings(as.numeric(as.character(.data[[pos_col]]))),
    requested_variant_id = .env$requested_ids, variant_id = NA_character_
  ) %>% .sv_target_output()

  loci <- targets %>% dplyr::group_by(query_chrom, query_pos) %>% dplyr::group_split()
  result <- purrr::map(loci, function(out) {
    tryCatch({
      chrom <- out$query_chrom[[1L]]
      pos <- out$query_pos[[1L]]
      if (!is.finite(pos) || pos < 1 || pos != floor(pos)) stop("Expected a positive integer VCF position.")
      region <- make_sv_region(chrom, pos, flank = 0L)
      out$region <- rep(region, nrow(out))
      raw <- extract_sv_region(vcf_file, region, tabix_path, allow_empty = TRUE) %>%
        dplyr::filter(CHROM == .env$chrom, as.numeric(POS) == .env$pos)

      purrr::map(seq_len(nrow(out)), function(i) {
        row <- out[i, ]
        tryCatch({
          id <- row$requested_variant_id[[1L]]
          record <- if (is.na(id)) raw else raw[which(raw$ID == id), , drop = FALSE]
          if (!nrow(record)) {
            row$status <- "not_found"
            row$error_message <- "No biallelic VCF record matches this position and requested ID."
            return(row)
          }
          if (nrow(record) > 1L) {
            row$status <- "ambiguous"
            row$error_message <- paste0("Multiple VCF records match. IDs: ", paste(unique(record$ID), collapse = ", "))
            return(row)
          }
          matched_id <- record$ID[[1L]]
          if (!is.na(matched_id) && !matched_id %in% c("", ".")) row$variant_id <- matched_id
          model_data <- .sv_target_model(
            phenotype_df, row$phenotype_id[[1L]], id_col, record, covariates, covariate_cols
          )
          .sv_fit_target_record(row, record, model_data, max_alt_af, max_missing, missing_as_reference)
        }, error = function(e) {
          row$error_message <- conditionMessage(e)
          row
        })
      }) %>% purrr::list_rbind()
    }, error = function(e) {
      out$error_message <- rep(conditionMessage(e), nrow(out))
      out
    })
  }, .progress = "SV sheet tests") %>% purrr::list_rbind() %>% arrange(sheet_row)
  result <- annotate_sv_results(result, phenotype_df, id_col)
  if (all(c("gene_id", "gene_symbol") %in% names(sv_sheet))) {
    result$gene_id <- coalesce(result$gene_id, as.character(sv_sheet$gene_id[result$sheet_row]))
    result$gene_symbol <- coalesce(result$gene_symbol, as.character(sv_sheet$gene_symbol[result$sheet_row]))
    result <- result %>% relocate(gene_id, gene_symbol, .after = phenotype_id)
  }
  if ("source_sheet_row" %in% names(sv_sheet)) {
    result$source_sheet_row <- sv_sheet$source_sheet_row[result$sheet_row]
  }
  result
}

# Strip Ensembl versions while preserving distinct PAR_Y identifiers.
.sv_stable_gene_id <- function(x) {
  sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+(?=(_PAR_Y)?$)", "\\1", trimws(as.character(x)), perl = TRUE)
}

.sv_gene_label <- function(x) {
  x <- trimws(as.character(x))
  x[!is.na(x) & x %in% c("", ".")] <- NA_character_
  x
}

# Normalize identifiers before matching. Keep unmapped inputs in mapping_qc.
harmonize_sv_inputs <- function(tss_locations, sv_sheet = NULL, expression = NULL,
                                proteomics = NULL, cds_gtf = NULL, sheet_col = NULL,
                                expression_col = "gene_id", proteomic_col = "Assay",
                                assay_map = NULL) {
  if (is.character(tss_locations) && length(tss_locations) == 1L) {
    tss_locations <- data.table::fread(file = tss_locations, check.names = FALSE)
  }
  tss <- tibble::as_tibble(tss_locations)
  if (!nrow(tss)) stop("The gene reference has no rows.")
  symbol_col <- intersect(c("gene_symbol", "gene_name"), names(tss))
  if (!"gene_id" %in% names(tss) || !length(symbol_col)) {
    stop("The gene reference must have gene_id and gene_name or gene_symbol columns.")
  }
  tss$gene_id <- .sv_stable_gene_id(tss$gene_id)
  if (anyNA(tss$gene_id) || any(!grepl("^ENS[A-Z]*G[0-9]+(_PAR_Y)?$", tss$gene_id))) {
    stop("The reference gene_id column must contain Ensembl gene IDs.")
  }
  tss$gene_name <- .sv_gene_label(tss[[symbol_col[[1L]]]])
  reference <- tss |> dplyr::distinct(gene_id, gene_symbol = gene_name)
  conflicts <- reference |> dplyr::filter(!is.na(gene_symbol)) |>
    dplyr::count(gene_id) |> dplyr::filter(n > 1L)
  if (nrow(conflicts)) {
    stop("Conflicting symbols in the reference for: ", paste(head(conflicts$gene_id, 5), collapse = ", "))
  }
  gene_map <- reference |> dplyr::arrange(is.na(gene_symbol)) |> dplyr::distinct(gene_id, .keep_all = TRUE)
  gene_index <- list2env(setNames(as.list(gene_map$gene_id), gene_map$gene_id), parent = emptyenv())
  symbol_index <- list2env(split(gene_map$gene_id, toupper(gene_map$gene_symbol)), parent = emptyenv())
  alias_index <- list()
  if (!is.null(assay_map)) {
    if (!all(c("Assay", "gene_id") %in% names(assay_map))) stop("assay_map needs Assay and gene_id columns.")
    alias_ids <- .sv_stable_gene_id(assay_map$gene_id)
    aliases <- .sv_gene_label(assay_map$Assay)
    if (anyNA(aliases) || any(!alias_ids %in% gene_map$gene_id)) {
      stop("Each assay_map row needs an assay label and an Ensembl ID in the reference.")
    }
    alias_index <- split(alias_ids, aliases)
  }

  resolve <- function(label, use_alias = FALSE) {
    if (is.na(label)) return(character())
    if (use_alias && length(alias_index[[label]])) return(unique(alias_index[[label]]))
    id <- .sv_stable_gene_id(label)
    if (!is.null(gene_index[[id]])) return(id)
    symbol_index[[toupper(label)]]
  }

  prepare <- function(data, id_col, dataset) {
    if (is.null(data)) return(list(data = NULL, qc = NULL))
    data <- tibble::as_tibble(data)
    if (length(id_col) != 1L || is.na(id_col) || !id_col %in% names(data)) {
      stop("Specify an existing identifier column for ", dataset, ".")
    }
    primary <- .sv_gene_label(data[[id_col]])
    other_cols <- intersect(setdiff(c("gene_id", "gene_symbol", "gene_name"), id_col), names(data))
    other <- lapply(data[other_cols], .sv_gene_label)
    mapped <- purrr::map(seq_len(nrow(data)), function(i) {
      extra <- vapply(other, function(x) x[[i]], character(1))
      labels <- c(primary[[i]], unname(extra))
      candidates <- c(list(resolve(labels[[1L]], dataset == "proteomics")), lapply(extra, resolve))
      unresolved_id <- !lengths(candidates) & grepl("^ENS[A-Z]*G[0-9]+", labels)
      known <- candidates[lengths(candidates) > 0L]
      if (!length(known) || any(unresolved_id)) return(c(gene_id = NA_character_, status = "unmapped"))
      ids <- Reduce(intersect, known)
      if (!length(ids)) return(c(gene_id = NA_character_, status = "conflict"))
      if (length(ids) > 1L) return(c(gene_id = NA_character_, status = "ambiguous"))
      c(gene_id = ids, status = "ok")
    })
    ids <- purrr::map_chr(mapped, "gene_id")
    status <- purrr::map_chr(mapped, "status")
    if (dataset %in% c("expression", "proteomics")) {
      phenotype_ids <- if (dataset == "expression") ids else primary
      status[status == "ok" & is.na(phenotype_ids)] <- "missing_phenotype_id"
      matched <- which(status == "ok")
      keys <- phenotype_ids[matched]
      status[matched[duplicated(keys) | duplicated(keys, fromLast = TRUE)]] <- "duplicate_phenotype"
    }
    symbols <- gene_map$gene_symbol[match(ids, gene_map$gene_id)]
    keep <- status == "ok"
    out <- data[keep, , drop = FALSE]
    if (dataset == "expression" && id_col != "gene_id") out[[id_col]] <- NULL
    if (dataset == "proteomics") {
      if (!id_col %in% c("Assay", "gene_id")) out[[id_col]] <- NULL
      out$Assay <- primary[keep]
    }
    out$gene_id <- ids[keep]
    out$gene_symbol <- symbols[keep]
    if (dataset == "sv_sheet" && !"source_sheet_row" %in% names(out)) out$source_sheet_row <- which(keep)
    qc <- tibble::tibble(dataset = dataset, input_row = seq_len(nrow(data)), input_label = primary,
                        gene_id = ids, gene_symbol = symbols, status = status)
    list(data = out, qc = qc)
  }

  if (!is.null(sv_sheet) && is.null(sheet_col)) {
    sheet_col <- intersect(c("gene_id", "gene_symbol", "gene_name"), names(sv_sheet))[1L]
  }
  sheet <- prepare(sv_sheet, sheet_col, "sv_sheet")
  expr <- prepare(expression, expression_col, "expression")
  protein <- prepare(proteomics, proteomic_col, "proteomics")
  if (!is.null(cds_gtf)) {
    if (is.character(cds_gtf) && length(cds_gtf) == 1L) cds_gtf <- read_cds_gtf(cds_gtf)
    cds_gtf <- tibble::as_tibble(cds_gtf)
    if (!all(c("gene_id", "chrom", "start", "end") %in% names(cds_gtf))) {
      stop("Supply a CDS table with gene_id, chrom, start, and end columns.")
    }
    if (!"gtf_gene_id" %in% names(cds_gtf)) cds_gtf$gtf_gene_id <- as.character(cds_gtf$gene_id)
    cds_gtf$gene_id <- .sv_stable_gene_id(cds_gtf$gtf_gene_id)
    cds_gtf$gene_symbol <- gene_map$gene_symbol[match(cds_gtf$gene_id, gene_map$gene_id)]
    mixed <- cds_gtf |> dplyr::distinct(gene_id, gtf_gene_id) |> dplyr::count(gene_id) |> dplyr::filter(n > 1L)
    if (nrow(mixed)) stop("The CDS table contains multiple versions of the same stable gene ID.")
  }
  list(
    gene_map = gene_map, tss_locations = tss, cds_gtf = cds_gtf,
    sv_sheet = sheet$data, expression = expr$data, proteomics = protein$data,
    mapping_qc = dplyr::bind_rows(sheet$qc, expr$qc, protein$qc)
  )
}
