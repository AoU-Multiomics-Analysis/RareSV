source("R/sv_set_tests.R")
source("R/sv_pgen_workflow.R")

run_smoke <- function() {
  fixtures <- tempfile("sv-skato-")
  dir.create(fixtures)
  on.exit(unlink(fixtures, recursive = TRUE), add = TRUE)
  sample_ids <- paste0("s", 1:25)
  gene_ids <- c("ENSG00000000001", "ENSG00000000002")
  custom_masks <- .sv_pgen_parse_masks("INS,DEL+INS")
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
    writer <- pgenlibr::NewPgenWriter(pgen, sample_ct = length(sample_ids), variant_ct = 6L)
    for (j in 1:6) {
      calls <- rep(if (j == 5L) -9L else 0L, length(sample_ids))
      calls[c(j, j + 7L)] <- 1L
      pgenlibr::AppendBiallelic(writer, calls)
    }
    pgenlibr::ClosePgenWriter(writer)
    writeLines(c("##fileformat=VCFv4.2",
                 '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
                 "#CHROM\tPOS\tID\tREF\tALT\tINFO",
                 sprintf("%d\t%d\tsv%d\tN\t<%s>\tEND=%d", chr, c(98:101, 147), 1:5,
                         c("DEL", "DEL", "DUP", "DUP", "DEL"), c(115:118, 160)),
                 sprintf("%d\t101\tsv6\tN\t<INS:ME:ALU>\t.", chr)), pvar)
    write_tsv(tibble(`#IID` = sample_ids), psam)
    scan <- function(...) {
      .sv_pgen_scan(pgen_path = pgen, pvar_path = pvar, psam_path = psam,
                    phenotype_bed_path = phenotype_bed, lof_path = lof,
                    chromosome = as.character(chr), max_alt_af = 0.1, ...)
    }
    list(exact = scan(window_padding_bp = 0L), padded = scan(window_padding_bp = 1L),
         default = scan(), reference = scan(missing_pgen_as_reference = TRUE),
         cds = scan(gtf_path = gtf),
         custom = scan(window_padding_bp = 0L, masks = custom_masks),
         custom_follow_up = scan(
           window_padding_bp = 0L, masks = custom_masks,
           significant_genes = tibble(gene_id = gene_ids[[1L]], chromosome = "chr1",
                                      p_value = 0.01, q_value = 0.02)
         ),
         follow_up = scan(
           missing_pgen_as_reference = TRUE,
           significant_genes = tibble(gene_id = gene_ids[[1L]], chromosome = "chr1",
                                      p_value = 0.01, q_value = 0.02)
         ))
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
  reference_masks <- map(scans, "reference") %>% list_rbind()
  stopifnot(
    identical(default_masks$n_svs_in_mask, rep(c(3L, 2L, 5L), 2L)),
    identical(default_masks$n_svs_eligible, rep(c(2L, 2L, 4L), 2L)),
    identical(reference_masks$n_svs_eligible, reference_masks$n_svs_in_mask),
    all(reference_masks$status == "ok"), all(is.finite(reference_masks$p_value))
  )
  masks <- map(scans, "cds") %>% list_rbind()
  genes <- .sv_pgen_gene_results(masks, "BH")
  stopifnot(
    nrow(masks) == 6L, all(masks$status == "ok"), all(masks$n_samples == 24L),
    all(masks$n_lof_carriers == 1L), all(masks$lof_covariate_used),
    identical(masks$n_svs_tested, rep(c(2L, 2L, 4L), 2)),
    nrow(genes) == 2L, all(is.finite(genes$p_value)), all(is.finite(genes$q_value))
  )
  follow_up <- map(scans, "follow_up") %>% list_rbind()
  stopifnot(
    nrow(follow_up) == 5L, all(follow_up$gene_id == gene_ids[[1L]]),
    setequal(follow_up$variant_id, paste0("sv", 1:5)),
    follow_up$tested_masks[follow_up$variant_id == "sv5"] == "DEL,DEL+DUP",
    all(follow_up$gene_q_value == 0.02), all(follow_up$n_samples == 24L),
    all(follow_up$lof_covariate_used), all(follow_up$status == "ok"),
    all(is.finite(follow_up$beta)), all(is.finite(follow_up$p_value)),
    nrow(scans[[2L]]$follow_up) == 0L
  )
  custom <- map(scans, "custom") %>% list_rbind()
  custom_genes <- .sv_pgen_gene_results(custom, "BH", masks = custom_masks)
  custom_follow_up <- map(scans, "custom_follow_up") %>% list_rbind()
  stopifnot(
    identical(custom$sv_mask, rep(c("INS", "DEL+INS"), 2L)),
    identical(custom$n_svs_in_mask, c(1L, 3L, 0L, 0L)),
    all(custom$n_svs_tested[1:2] == c(1L, 3L)),
    identical(custom_genes$status, c("ok", "no_test")),
    isTRUE(all.equal(custom_genes$p_value[[1L]], .sv_cauchy_p(custom$p_value[1:2]))),
    setequal(custom_follow_up$variant_id, c("sv1", "sv2", "sv6")),
    custom_follow_up$tested_masks[custom_follow_up$variant_id == "sv6"] == "INS,DEL+INS",
    all(custom_follow_up$status == "ok")
  )
  cat("PASS: split BEDs, BED and CDS overlaps, pgenlibr hard calls, and two-chromosome gathering.\n")
}
run_smoke()
