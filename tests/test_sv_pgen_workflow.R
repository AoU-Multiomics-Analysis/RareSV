source("R/sv_set_tests.R")
source("R/sv_pgen_workflow.R")

pvar <- tibble(
  CHROM = c("1", "1", "1"), POS = c(99L, 105L, 110L),
  ID = c("del", "dup", "padding_only"), REF = "N",
  ALT = c("<DEL>", "<DUP>", "<DEL>"), INFO = c("END=104", "END=120", "END=140")
)
cds <- tibble(gene_id = "g1", chrom = "chr1", start = c(103L, 109L), end = c(104L, 110L))
links <- .sv_pgen_interval_links(pvar, cds)
stopifnot(identical(links$pgen_index, c(1L, 2L)))

samples <- paste0("s", 1:7)
expression <- tibble(sample_id = samples, phenotype = c(1, 2, 0, 3, -1, 4, 5))
Z <- cbind(del = c(1L, 0L, 0L, 0L, 0L, 0L, 0L), dup = c(0L, 1L, NA, 0L, 0L, 0L, 0L))
rownames(Z) <- samples
lof <- tibble(sample_id = c("s1", "s7", "outside"), has_lof_variant = c(1L, NA_integer_, 1L))
test_env <- new.env()
sys.source("R/sv_pgen_workflow.R", envir = test_env)
null_calls <- 0L
mask_calls <- list()
test_env$.sv_fit_null_model <- function(model_data) {
  null_calls <<- null_calls + 1L
  stopifnot(identical(model_data$data$sample_id, samples[1:6]),
            identical(model_data$data$has_lof_variant, c(1L, rep(0L, 5))),
            identical(all.vars(model_data$formula), c("phenotype", "has_lof_variant")))
  "adjusted null"
}
test_env$.sv_run_skat <- function(Z, null_model, test, max_missing) {
  stopifnot(null_model == "adjusted null", test == "SKATO",
            identical(rownames(Z), samples[1:6]))
  mask_calls[[length(mask_calls) + 1L]] <<- colnames(Z)
  list(p.value = 0.02, param = list(n.marker.test = ncol(Z)))
}
result <- test_env$.sv_pgen_test_gene(
  expression, Z, c("DEL", "DUP"), lof, TRUE, NULL, character(), 0.2, 0.2
)
stopifnot(
  null_calls == 1L, identical(mask_calls, list("del", "dup", c("del", "dup"))),
  identical(result$sv_mask, c("DEL", "DUP", "DEL+DUP")),
  identical(result$n_svs_tested, c(1L, 1L, 2L)),
  all(result$n_samples == 6L), all(result$n_lof_carriers == 1L)
)
stopifnot(
  isTRUE(all.equal(.sv_cauchy_p(c(.01, .04)), .01600759, tolerance = 1e-6)),
  identical(.sv_cauchy_p(.03), .03),
  isTRUE(all.equal(.sv_cauchy_p(c(1e-100, .5)), 2e-100, tolerance = 1e-12, scale = 2e-100))
)

mask_rows <- tibble(
  gene_id = rep(paste0("g", 1:3), each = 3), gene_symbol = rep(LETTERS[1:3], each = 3),
  chromosome = rep(as.character(1:3), each = 3),
  sv_mask = rep(c("DEL", "DUP", "DEL+DUP"), 3),
  status = c("ok", "ok", "ok", "ok", "no_eligible_sv", "no_eligible_sv",
             "ok", "error", "ok"),
  p_value = c(.01, .01, .01, .04, NA, NA, .001, NA, .001)
)
gathered <- .sv_pgen_gene_results(mask_rows, "BH")
stopifnot(
  isTRUE(all.equal(gathered$p_value, c(.01, .04, NA))),
  isTRUE(all.equal(gathered$q_value, c(.03, .06, NA))),
  identical(gathered$status, c("ok", "ok", "incomplete")),
  identical(gathered$significant, c(TRUE, FALSE, FALSE))
)
stopifnot(inherits(try(.sv_pgen_gene_results(bind_rows(mask_rows, mask_rows), "BH"), silent = TRUE), "try-error"))

uniform_p <- seq(0.005, 0.995, length.out = 100)
storey_rows <- tibble(
  gene_id = rep(paste0("uniform_", seq_along(uniform_p)), each = 3),
  gene_symbol = NA_character_, chromosome = rep(c("1", "2"), each = 150),
  sv_mask = rep(names(.sv_pgen_masks), length(uniform_p)),
  status = "ok", p_value = rep(uniform_p, each = 3)
)
storey <- .sv_pgen_gene_results(storey_rows)
stopifnot(all(storey$fdr_method == "storey"), all(is.finite(storey$q_value)),
          isTRUE(all.equal(storey$q_value, qvalue::qvalue(uniform_p, lfdr.out = FALSE)$qvalues)))

test_pgen_cli <- function() {
  script <- normalizePath("R/sv_pgen_workflow.R")
  runner <- file.path(R.home("bin"), "Rscript")
  work <- tempfile("sv-skato-cli-")
  dir.create(work)
  previous <- setwd(work)
  on.exit({setwd(previous); unlink(work, recursive = TRUE)}, add = TRUE)
  rows <- bind_cols(
    dplyr::select(mask_rows[1:3, ], gene_id, gene_symbol, chromosome),
    .sv_pgen_mask_rows() %>% mutate(status = "ok", p_value = 0.01)
  )
  write_tsv(rows, "first shard's.tsv")
  write_tsv(mutate(rows, gene_id = "g3", status = c("ok", "error", "ok"),
                   p_value = c(0.001, NA, 0.001),
                   error_message = c(NA, "synthetic failure", NA)), "second.tsv")
  writeLines(file.path(work, c("first shard's.tsv", "second.tsv")), "files.txt")
  run <- function(...) {
    suppressWarnings(system2(runner, shQuote(c(script, ...)), stdout = TRUE, stderr = TRUE))
  }
  result <- run("gather", "--mask-list", "files.txt", "--fdr-method", "BH")
  stopifnot(is.null(attr(result, "status")))
  genes <- read_tsv("gene_results.tsv", show_col_types = FALSE)
  stopifnot(identical(genes$q_value, c(0.02, NA_real_)))
  unresolved <- run("scan", "--pgen", "gs://test/genotypes.pgen")
  stopifnot(attr(unresolved, "status") != 0L,
            any(grepl("Localization error for PGEN: unresolved URI", unresolved, fixed = TRUE)))
  writeLines("", "input")
  unresolved <- run("scan", "--pgen", "input", "--pvar", "input", "--psam", "input",
                    "--phenotype-bed", "gs://test/expression.bed.gz", "--lof", "input")
  stopifnot(attr(unresolved, "status") != 0L,
            any(grepl("Localization error for phenotype BED: unresolved URI", unresolved, fixed = TRUE)))
  unresolved <- run("scan", "--pgen", "input", "--pvar", "input", "--psam", "input",
                    "--phenotype-bed", "input", "--lof", "input", "--gtf", "gs://test/cds.gtf")
  stopifnot(attr(unresolved, "status") != 0L,
            any(grepl("Localization error for CDS GTF: unresolved URI", unresolved, fixed = TRUE)))
  writeLines("gs://test/mask.tsv", "files.txt")
  unresolved <- run("gather", "--mask-list", "files.txt")
  stopifnot(attr(unresolved, "status") != 0L,
            any(grepl("Localization error for mask results: unresolved URI", unresolved, fixed = TRUE)))
}
test_pgen_cli()

test_phenotype_bed <- function() {
  path <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(path))
  bed <- tibble(`#chr` = "chr1", start = 0, end = 1, gene_id = "ENSG00000000001.2",
                `000100` = -1.5, `1000291` = NA_real_)
  write_tsv(bed, path)
  parsed <- .sv_pgen_read_phenotype_bed(path, "1")
  stopifnot(identical(names(parsed), names(bed)), parsed$gene_id == "ENSG00000000001.2",
            parsed$start == 0, parsed$end == 1, parsed[["000100"]] == -1.5,
            is.na(parsed[["1000291"]]))

  expect_bed_error <- function(input, message) {
    write_tsv(input, path)
    result <- tryCatch(.sv_pgen_read_phenotype_bed(path, "1"), error = conditionMessage)
    stopifnot(is.character(result), grepl(message, result, fixed = TRUE))
  }
  expect_bed_error(select(bed, gene_id, `000100`), "must begin with")
  expect_bed_error(mutate(bed, `#chr` = "chr2"), "only rows on chromosome")
  expect_bed_error(mutate(bed, start = -1), "0 <= start < end")
  expect_bed_error(bind_rows(bed, mutate(bed, gene_id = "ENSG00000000001.3")),
                   "one row per stable gene ID")
}
test_phenotype_bed()
cat("PASS: CDS overlap, LoF-adjusted mask tests, Cauchy combination, and global FDR.\n")
