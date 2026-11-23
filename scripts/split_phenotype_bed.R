#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/split_phenotype_bed.R <expression.tsv> <tss.tsv> <output_dir>",
       call. = FALSE)
}

# Some expression exports have a metadata line above the sample headers.
expression <- data.table::fread(file = args[[1L]], skip = "gene_id", header = TRUE,
                                check.names = FALSE, showProgress = FALSE) |>
  as_tibble()
tss <- read_tsv(args[[2L]], show_col_types = FALSE)
out_dir <- args[[3L]]

if (!"gene_id" %in% names(expression) || ncol(expression) < 2L) {
  stop("The expression table needs gene_id and sample columns.")
}
if (!all(c("gene_id", "chrom", "tss_1based") %in% names(tss))) {
  stop("The TSS table needs gene_id, chrom, and tss_1based columns.")
}

stable_id <- function(x) {
  sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+(?=(_PAR_Y)?$)", "\\1", trimws(x), perl = TRUE)
}

locations <- tss |>
  transmute(.gene_key = stable_id(gene_id), .chrom = trimws(as.character(chrom)), .tss = tss_1based) |>
  filter(!is.na(.chrom), .chrom != "") |>
  distinct()
mapped <- expression |>
  mutate(.gene_key = stable_id(gene_id)) |>
  left_join(locations, by = ".gene_key", relationship = "many-to-one", na_matches = "never")

unmatched <- sum(is.na(mapped$.chrom))
if (unmatched == nrow(mapped)) stop("No expression rows matched the TSS table.")
if (unmatched > 0L) warning("Skipping ", unmatched, " rows with no chromosome match.", call. = FALSE)
mapped <- filter(mapped, !is.na(.chrom))
if (anyDuplicated(mapped$.gene_key)) stop("Use one expression row per stable gene ID.")
if (any(!str_detect(mapped$.chrom, "^[[:alnum:]][[:alnum:]_.-]*$"))) {
  stop("A matched gene has an invalid chromosome name.")
}
if (!is.numeric(mapped$.tss) ||
    any(!is.finite(mapped$.tss) | mapped$.tss < 1 | mapped$.tss != floor(mapped$.tss))) {
  stop("A matched gene has an invalid TSS. Use positive integer positions.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
value_cols <- setdiff(names(expression), c("gene_id", "gene_name"))
mapped |>
  group_by(.chrom) |>
  group_walk(function(rows, key) {
    chr <- key$.chrom
    bed <- rows |>
      arrange(.tss, gene_id) |>
      transmute("#chr" = chr, start = .tss - 1L, end = .tss, gene_id, across(all_of(value_cols)))
    path <- file.path(out_dir, paste0("expression_", chr, ".bed.gz"))
    write_tsv(bed, path)
    message("Wrote ", nrow(bed), " genes to ", path)
  })
