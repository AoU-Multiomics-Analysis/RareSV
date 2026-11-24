source("R/sv_set_tests.R")
source("R/sv_pgen_workflow.R")

run_smoke <- function() {
  fixtures <- tempfile("sv-skato-")
  dir.create(fixtures)
  on.exit(unlink(fixtures, recursive = TRUE), add = TRUE)
  sample_ids <- paste0("s", 1:25)
  gene_ids <- c("ENSG00000000001", "ENSG00000000002")
  gtf <- file.path(fixtures, "coding.gtf")
  writeLines(sprintf('chr%d\tsmoke\tCDS\t100\t140\t.\t+\t0\tgene_id "%s.1";',
                     1:2, gene_ids), gtf)
  lof <- file.path(fixtures, "lof.tsv.gz")
  write_tsv(tibble(sample_id = c("s6", "s13", "s25"), gene_id = paste0(gene_ids[c(1, 2, 2)], ".1"),
                   has_lof_variant = TRUE), lof)
  scores <- as_tibble(set_names(as.list(sin(1:24)), rev(sample_ids[1:24])))
  expression_path <- file.path(fixtures, "expression scores.tsv")
  writeLines("# Expression score export", expression_path)
  expression <- bind_cols(tibble(gene_id = gene_ids, phenotype_only = 0), scores[rep(1L, 2L), ])
  cat(paste(format_tsv(expression), collapse = ""), file = expression_path, append = TRUE)
  tss_path <- file.path(fixtures, "tss.tsv")
  write_tsv(tibble(gene_id = paste0(gene_ids, ".1"), chrom = paste0("chr", 1:2),
                   tss_1based = c(101L, 400L)), tss_path)
  splitter <- c("scripts/split_phenotype_bed.R", expression_path, tss_path, fixtures)
  stopifnot(system2(file.path(R.home("bin"), "Rscript"), shQuote(splitter)) == 0L)
  scans <- map(seq_along(gene_ids), function(chr) {
    pgen <- file.path(fixtures, paste0(chr, ".pgen"))
    pvar <- file.path(fixtures, paste0(chr, ".pvar"))
    psam <- file.path(fixtures, paste0(chr, ".psam"))
    phenotype_bed <- file.path(fixtures, paste0("expression_chr", chr, ".bed.gz"))
    bed <- read_tsv(phenotype_bed, show_col_types = FALSE)
    stopifnot(bed$start == c(100L, 399L)[chr], bed$end == bed$start + 1L)
    writer <- pgenlibr::NewPgenWriter(pgen, sample_ct = length(sample_ids), variant_ct = 5L)
    for (j in 1:5) {
      calls <- integer(length(sample_ids))
      calls[c(j, j + 7L)] <- 1L
      pgenlibr::AppendBiallelic(writer, calls)
    }
    pgenlibr::ClosePgenWriter(writer)
    writeLines(c("##fileformat=VCFv4.2",
                 '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
                 "#CHROM\tPOS\tID\tREF\tALT\tINFO",
                 sprintf("%d\t%d\tsv%d\tN\t<%s>\tEND=%d", chr, c(98:101, 147), 1:5,
                         c("DEL", "DEL", "DUP", "DUP", "DEL"), c(115:118, 160))), pvar)
    write_tsv(tibble(`#IID` = sample_ids), psam)
    scan <- function(...) {
      .sv_pgen_scan(pgen_path = pgen, pvar_path = pvar, psam_path = psam,
                    phenotype_bed_path = phenotype_bed, lof_path = lof,
                    chromosome = as.character(chr), max_alt_af = 0.1, ...)
    }
    list(exact = scan(window_padding_bp = 0L), padded = scan(window_padding_bp = 1L),
         default = scan(), cds = scan(gtf_path = gtf))
  })
  bed_masks <- map(scans, "exact") %>% list_rbind()
  bed_genes <- .sv_pgen_gene_results(bed_masks, "BH")
  stopifnot(
    identical(bed_masks$n_svs_in_mask, c(2L, 1L, 3L, 0L, 0L, 0L)),
    all(bed_masks$status[4:6] == "no_overlapping_sv"),
    identical(bed_genes$status, c("ok", "no_test"))
  )
  padded_masks <- map(scans, "padded") %>% list_rbind()
  stopifnot(identical(padded_masks$n_svs_in_mask, c(2L, 2L, 4L, 0L, 0L, 0L)))
  default_masks <- map(scans, "default") %>% list_rbind()
  stopifnot(identical(default_masks$n_svs_in_mask, rep(c(3L, 2L, 5L), 2L)))
  masks <- map(scans, "cds") %>% list_rbind()
  genes <- .sv_pgen_gene_results(masks, "BH")
  stopifnot(
    nrow(masks) == 6L, all(masks$status == "ok"), all(masks$n_samples == 24L),
    all(masks$n_lof_carriers == 1L), all(masks$lof_covariate_used),
    identical(masks$n_svs_tested, rep(c(2L, 2L, 4L), 2)),
    nrow(genes) == 2L, all(is.finite(genes$p_value)), all(is.finite(genes$q_value))
  )
  cat("PASS: split BEDs, BED and CDS overlaps, pgenlibr hard calls, and two-chromosome gathering.\n")
}
run_smoke()
