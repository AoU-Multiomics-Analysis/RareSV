# Source sv_set_tests.R before using these functions from an R session.
.sv_pgen_masks <- list(DEL = "DEL", DUP = "DUP", `DEL+DUP` = c("DEL", "DUP"))

.sv_pgen_chr <- function(chromosome) {
  x <- str_remove(str_to_upper(trimws(as.character(chromosome))), "^CHR")
  recode(x, "23" = "X", "24" = "Y", "26" = "MT", "M" = "MT")
}

.sv_local_input <- function(path, label) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Supply a path for ", label, ".")
  }
  if (str_detect(path, "^[[:alpha:]][[:alnum:]+.-]*://")) {
    stop("Localization error for ", label, ": unresolved URI ", path)
  }
  if (!file.exists(path) || dir.exists(path) || file.access(path, 4L) != 0L) {
    stop("Localization error for ", label, ": file is not readable: ", path)
  }
  path
}

.sv_pgen_cds_links <- function(pvar, cds) {
  spans <- .sv_cds_spans(pvar) %>%
    mutate(pgen_index = row_number(), chrom = .sv_pgen_chr(CHROM)) %>%
    filter(sv_type %in% c("DEL", "DUP"), !is.na(ALT), !str_detect(ALT, ","))
  unresolved <- sum(spans$span_status == "unresolved_span")
  if (unresolved) {
    if (!any(spans$span_status == "ok")) {
      stop("No DEL/DUP spans could be resolved from the PVAR. Preserve INFO/END for symbolic SVs.")
    }
    warning("Excluded ", unresolved, " DEL/DUP records with unresolved spans. Check PVAR INFO/END.",
            call. = FALSE)
  }
  spans <- filter(spans, span_status == "ok")
  cds <- mutate(cds, chrom = .sv_pgen_chr(chrom))
  links <- map(intersect(unique(cds$chrom), unique(spans$chrom)), function(contig) {
    coding <- filter(cds, chrom == .env$contig)
    variants <- filter(spans, chrom == .env$contig)
    hits <- IRanges::findOverlaps(
      IRanges::IRanges(coding$start, coding$end),
      IRanges::IRanges(variants$affected_start, variants$affected_end)
    )
    tibble(gene_id = coding$gene_id[S4Vectors::queryHits(hits)],
           pgen_index = variants$pgen_index[S4Vectors::subjectHits(hits)],
           sv_type = variants$sv_type[S4Vectors::subjectHits(hits)])
  })
  bind_rows(tibble(gene_id = character(), pgen_index = integer(), sv_type = character()), links) %>%
    distinct() %>% arrange(gene_id, pgen_index)
}

.sv_pgen_mask_rows <- function() {
  tibble(sv_mask = names(.sv_pgen_masks), status = "no_cds_sv", error_message = NA_character_,
         n_samples = NA_integer_, n_lof_missing = NA_integer_, n_lof_carriers = NA_integer_,
         lof_covariate_used = NA, n_svs_in_mask = 0L, n_svs_eligible = 0L,
         n_svs_tested = 0L, p_value = NA_real_, warnings = NA_character_)
}

.sv_pgen_test_gene <- function(pheno, Z, sv_types, lof_status, lof_missing_as_noncarrier,
                               covariates, covariate_cols, max_alt_af, max_missing) {
  rows <- .sv_pgen_mask_rows()
  if (!ncol(Z)) return(rows)
  model <- .sv_model_data(pheno, covariates, covariate_cols)
  model <- .sv_with_lof_covariate(model, lof_status, lof_missing_as_noncarrier)
  analysis_ids <- model$data$sample_id
  if (length(analysis_ids) < 3L) stop("Fewer than three samples remain after covariate filtering.")
  stopifnot(all(analysis_ids %in% rownames(Z)), length(sv_types) == ncol(Z))
  Z <- Z[analysis_ids, , drop = FALSE]
  filtered <- filter_sv_genotypes(Z, max_alt_af, max_missing)
  keep <- filtered$variant_qc$keep
  null_model <- if (any(keep)) .sv_fit_null_model(model)
  rows <- rows %>%
    mutate(n_samples = length(.env$analysis_ids),
           n_lof_missing = model$lof_qc$n_lof_missing,
           n_lof_carriers = model$lof_qc$n_lof_carriers,
           lof_covariate_used = model$lof_qc$lof_covariate_used)
  map(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    in_mask <- sv_types %in% .sv_pgen_masks[[row$sv_mask]]
    eligible <- in_mask & keep
    row$n_svs_in_mask <- sum(in_mask)
    row$n_svs_eligible <- sum(eligible)
    if (!any(in_mask)) return(row)
    row$status <- "no_eligible_sv"
    if (!any(eligible)) return(row)
    warnings <- character()
    row <- tryCatch(withCallingHandlers({
      fit <- .sv_run_skat(Z[, eligible, drop = FALSE], null_model, "SKATO", max_missing)
      p <- fit$p.value
      n_tested <- fit$param$n.marker.test
      if (length(n_tested) != 1L || !is.finite(n_tested) || n_tested < 1L ||
          length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
        stop("SKAT-O returned an invalid p-value or tested zero SVs.")
      }
      mutate(row, status = "ok", n_svs_tested = as.integer(n_tested), p_value = p)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) mutate(row, status = "error", error_message = conditionMessage(e)))
    if (length(warnings)) row$warnings <- str_c(unique(warnings), collapse = " | ")
    row
  }) %>% list_rbind()
}

.sv_cauchy_p <- function(p) {
  if (!length(p)) return(NA_real_)
  stopifnot(is.numeric(p), all(is.finite(p) & p >= 0 & p <= 1))
  if (length(p) == 1L || all(p == p[[1L]])) return(p[[1L]])
  # Keep finite Cauchy quantiles at the two numerical boundaries.
  bounded <- pmin(pmax(p, .Machine$double.xmin), 1 - .Machine$double.eps)
  stats::pcauchy(mean(stats::qcauchy(bounded, lower.tail = FALSE)), lower.tail = FALSE)
}

.sv_pgen_gene_results <- function(raw, fdr_method = c("storey", "BH")) {
  fdr_method <- match.arg(fdr_method)
  allowed <- c("ok", "no_cds_sv", "no_eligible_sv", "error")
  if (anyNA(raw$gene_id) || any(!raw$sv_mask %in% names(.sv_pgen_masks)) ||
      any(!raw$status %in% allowed) || anyDuplicated(raw[c("gene_id", "sv_mask")])) {
    stop("Invalid or duplicate gene-mask rows. Each gene must belong to one chromosome shard.")
  }
  if (any(!is.finite(raw$p_value[raw$status == "ok"]) |
          raw$p_value[raw$status == "ok"] < 0 | raw$p_value[raw$status == "ok"] > 1)) {
    stop("Successful mask rows must contain valid p-values.")
  }
  genes <- raw %>%
    group_by(gene_id, gene_symbol, chromosome) %>%
    summarise(
      n_masks = n(), n_masks_tested = sum(status == "ok"),
      n_masks_failed = sum(status == "error"),
      p_value = if (n_masks_failed == 0L) .sv_cauchy_p(p_value[status == "ok"]) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(status = case_when(n_masks_failed > 0L ~ "incomplete",
                              n_masks_tested == 0L ~ "no_test", TRUE ~ "ok"),
           q_value = NA_real_)
  if (any(genes$n_masks != length(.sv_pgen_masks)) || anyDuplicated(genes$gene_id)) {
    stop("Expected all three mask rows for each gene.")
  }
  family <- which(genes$status != "no_test")
  if (any(genes$status == "ok")) {
    # Keep incomplete tests in the family at p=1. Report their q-values as missing.
    p <- coalesce(genes$p_value[family], 1)
    adjusted <- if (fdr_method == "storey") .sv_storey_qvalues(p) else stats::p.adjust(p, method = "BH")
    genes$q_value[family] <- if_else(genes$status[family] == "ok", adjusted, NA_real_)
  }
  genes %>% mutate(fdr_method = .env$fdr_method, significant = coalesce(q_value < 0.05, FALSE)) %>%
    arrange(is.na(p_value), p_value, gene_id)
}

.sv_pgen_read_tsv <- function(path, ..., colClasses = "character", sep = "\t") {
  data.table::fread(file = path, sep = sep, header = TRUE, quote = "", check.names = FALSE,
                   colClasses = colClasses, showProgress = FALSE, ...) %>%
    as_tibble()
}

.sv_pgen_scan <- function(pgen_path, pvar_path, psam_path, trait_path, gtf_path, lof_path, chromosome,
                          covariate_path = NULL, covariate_cols = character(),
                          max_alt_af = 0.01, max_missing = 0.05, lof_missing_as_noncarrier = TRUE) {
  stopifnot(length(chromosome) == 1L, !is.na(chromosome), nzchar(chromosome),
            is.finite(max_alt_af), max_alt_af > 0, max_alt_af <= 1,
            is.finite(max_missing), max_missing >= 0, max_missing <= 1)
  chr <- .sv_pgen_chr(chromosome)
  pvar <- .sv_pgen_read_tsv(pvar_path, skip = "#CHROM")
  if (!all(c("#CHROM", "POS", "ID", "REF", "ALT") %in% names(pvar))) {
    stop("PVAR must have #CHROM, POS, ID, REF, and ALT columns.")
  }
  pvar <- rename(pvar, CHROM = `#CHROM`)
  if (!nrow(pvar) || anyNA(pvar$CHROM) || any(.sv_pgen_chr(pvar$CHROM) != chr)) {
    stop("PVAR must contain only variants on chromosome ", chromosome, ".")
  }
  if (!"INFO" %in% names(pvar)) pvar$INFO <- "."
  if (anyNA(pvar$POS) || any(!str_detect(pvar$POS, "^[1-9][0-9]*$")) ||
      any(as.double(pvar$POS) > .Machine$integer.max)) stop("PVAR contains invalid positions.")

  psam <- .sv_pgen_read_tsv(psam_path, sep = "auto")
  if ("#IID" %in% names(psam)) psam <- rename(psam, IID = `#IID`)
  if (!"IID" %in% names(psam) || !nrow(psam) || anyNA(.sv_gene_label(psam$IID)) ||
      anyDuplicated(psam$IID)) stop("PSAM must have unique, nonblank IID values.")
  traits <- .sv_pgen_read_tsv(trait_path, colClasses = list(character = "gene_id"))
  if (!"gene_id" %in% names(traits)) stop("The molecular trait file needs a gene_id column.")
  if (anyNA(.sv_gene_label(traits$gene_id))) stop("Molecular trait gene IDs must be nonblank.")
  sample_index <- which(psam$IID %in% names(traits))
  if (!length(sample_index)) stop("No PSAM IIDs match molecular trait sample columns.")
  sample_ids <- psam$IID[sample_index]
  trait_keys <- .sv_stable_gene_id(traits$gene_id)
  cds <- read_cds_gtf(gtf_path) %>%
    mutate(gene_id = .sv_stable_gene_id(gene_id)) %>%
    filter(gene_id %in% .env$trait_keys, .sv_pgen_chr(chrom) == .env$chr)
  selected <- which(trait_keys %in% cds$gene_id)
  if (!length(selected)) stop("No molecular traits match CDS genes on chromosome ", chromosome, ".")
  if (anyDuplicated(trait_keys[selected])) stop("Use one molecular trait row per stable gene ID.")
  symbols <- if ("gene_symbol" %in% names(traits)) .sv_gene_label(traits$gene_symbol) else
    rep(NA_character_, nrow(traits))
  if ("gene_symbol" %in% names(cds)) {
    cds_symbols <- cds %>% distinct(gene_id, gene_symbol) %>% filter(!is.na(gene_symbol))
    symbols <- coalesce(symbols, cds_symbols$gene_symbol[match(trait_keys, cds_symbols$gene_id)])
  }
  links <- .sv_pgen_cds_links(pvar, cds)
  if (anyNA(pvar$ID[links$pgen_index]) || any(pvar$ID[links$pgen_index] %in% c("", ".")) ||
      anyDuplicated(pvar$ID[unique(links$pgen_index)])) {
    stop("CDS-overlapping SVs must have unique, nonblank PVAR IDs.")
  }
  links_by_gene <- split(seq_len(nrow(links)), links$gene_id)
  lof <- .sv_prepare_lof_table(
    .sv_pgen_read_tsv(lof_path, select = c("sample_id", "gene_id", "has_lof_variant")),
    lof_missing_as_noncarrier
  ) %>%
    filter(gene_id %in% .env$trait_keys[.env$selected])
  lof_by_gene <- split(seq_len(nrow(lof)), lof$gene_id)
  covariates <- if (!is.null(covariate_path)) {
    .sv_pgen_read_tsv(covariate_path, colClasses = list(character = "sample_id"))
  }
  pvar_handle <- pgenlibr::NewPvar(pvar_path)
  on.exit(pgenlibr::ClosePvar(pvar_handle), add = TRUE)
  pgen <- pgenlibr::NewPgen(pgen_path, pvar = pvar_handle, sample_subset = sample_index)
  on.exit(pgenlibr::ClosePgen(pgen), add = TRUE, after = FALSE)
  if (pgenlibr::GetRawSampleCt(pgen) != nrow(psam) || pgenlibr::GetVariantCt(pgen) != nrow(pvar)) {
    stop("PGEN dimensions do not match the PSAM and PVAR.")
  }
  message("Chromosome ", chromosome, ": ", length(selected), " genes and ",
          length(sample_ids), " matched sample columns.")

  map(selected, function(index) {
    id <- trait_keys[[index]]
    candidate <- links[links_by_gene[[id]], , drop = FALSE]
    result <- .sv_pgen_mask_rows()
    if (nrow(candidate)) {
      result <- tryCatch({
        pheno <- suppressMessages(prep_burden_phenotype(
          traits[index, ], traits$gene_id[[index]], "gene_id", sample_ids
        ))
        Z <- pgenlibr::ReadIntList(pgen, candidate$pgen_index)
        dimnames(Z) <- list(sample_ids, pvar$ID[candidate$pgen_index])
        suppressMessages(.sv_pgen_test_gene(
          pheno, Z, candidate$sv_type, lof[lof_by_gene[[id]], , drop = FALSE],
          lof_missing_as_noncarrier, covariates, covariate_cols, max_alt_af, max_missing
        ))
      }, error = function(e) mutate(result, status = "error", error_message = conditionMessage(e)))
    }
    result %>% mutate(gene_id = .env$id, gene_symbol = .env$symbols[[index]],
                      chromosome = .env$chr, .before = 1L)
  }) %>% list_rbind()
}

.sv_pgen_main <- function(args) {
  option <- optparse::make_option
  if (!length(args) || !args[[1L]] %in% c("scan", "gather")) stop("Specify scan or gather.")
  if (args[[1L]] == "scan") {
    options <- list(
      option("--pgen"), option("--pvar"), option("--psam"), option("--traits"), option("--gtf"),
      option("--lof"), option("--chromosome"), option("--covariates"), option("--covariate-cols"),
      option("--max-alt-af", type = "double", default = 0.01),
      option("--max-missing", type = "double", default = 0.05),
      option("--unknown-unlisted-lof", action = "store_true", default = FALSE)
    )
    opt <- optparse::parse_args(optparse::OptionParser(option_list = options), args = args[-1L],
                               convert_hyphens_to_underscores = TRUE)
    paths <- imap_chr(list(PGEN = opt$pgen, PVAR = opt$pvar, PSAM = opt$psam,
                           traits = opt$traits, GTF = opt$gtf, LoF = opt$lof), .sv_local_input)
    cov_path <- if (!is.null(opt$covariates)) .sv_local_input(opt$covariates, "covariates")
    cov_cols <- if (!is.null(opt$covariate_cols)) str_trim(str_split_1(opt$covariate_cols, ",")) else character()
    if (xor(!is.null(cov_path), length(cov_cols) > 0L)) {
      stop("Supply both a covariate file and the covariate column names.")
    }
    results <- .sv_pgen_scan(
      paths[["PGEN"]], paths[["PVAR"]], paths[["PSAM"]], paths[["traits"]], paths[["GTF"]],
      paths[["LoF"]], opt$chromosome, cov_path, cov_cols, opt$max_alt_af, opt$max_missing,
      !opt$unknown_unlisted_lof
    )
    readr::write_tsv(results, "mask_results.tsv")
    message("Wrote ", nrow(results), " gene-mask rows.")
  } else {
    options <- list(option("--mask-list"), option("--fdr-method", default = "storey"))
    opt <- optparse::parse_args(optparse::OptionParser(option_list = options), args = args[-1L],
                               convert_hyphens_to_underscores = TRUE)
    paths <- readLines(.sv_local_input(opt$mask_list, "mask result file list"), warn = FALSE)
    if (!length(paths)) stop("The mask result file list is empty.")
    paths <- map_chr(paths, .sv_local_input, label = "mask results")
    raw <- map(paths, function(path) {
      readr::read_tsv(path, col_types = readr::cols(
        .default = readr::col_guess(), chromosome = readr::col_character(),
        gene_id = readr::col_character(), gene_symbol = readr::col_character(),
        error_message = readr::col_character(), warnings = readr::col_character(),
        p_value = readr::col_double()
      ), show_col_types = FALSE)
    }) %>% list_rbind()
    genes <- .sv_pgen_gene_results(raw, opt$fdr_method)
    readr::write_tsv(raw, "mask_results.tsv")
    readr::write_tsv(genes, "gene_results.tsv")
    readr::write_tsv(filter(genes, significant), "significant_genes.tsv")
    message("Gathered ", nrow(genes), " genes and found ", sum(genes$significant), " discoveries at 5% FDR.")
  }
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
  if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
  source(file.path(dirname(script), "sv_set_tests.R"))
  .sv_pgen_main(commandArgs(trailingOnly = TRUE))
}
