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
  lof <- file.path(fixtures, "lof.tsv")
  write_tsv(tibble(sample_id = c("s6", "s13", "s25"), gene_id = paste0(gene_ids[c(1, 2, 2)], ".1"),
                   has_lof_variant = TRUE), lof)
  masks <- map(seq_along(gene_ids), function(chr) {
    pgen <- file.path(fixtures, paste0(chr, ".pgen"))
    pvar <- file.path(fixtures, paste0(chr, ".pvar"))
    psam <- file.path(fixtures, paste0(chr, ".psam"))
    traits <- file.path(fixtures, paste0(chr, ".tsv"))
    writer <- pgenlibr::NewPgenWriter(pgen, sample_ct = length(sample_ids), variant_ct = 4L)
    for (j in 1:4) {
      calls <- integer(length(sample_ids))
      calls[c(j, j + 7L)] <- 1L
      pgenlibr::AppendBiallelic(writer, calls)
    }
    pgenlibr::ClosePgenWriter(writer)
    writeLines(c("##fileformat=VCFv4.2",
                 '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
                 "#CHROM\tPOS\tID\tREF\tALT\tINFO",
                 sprintf("%d\t%d\tsv%d\tN\t<%s>\tEND=%d", chr, 98:101, 1:4,
                         c("DEL", "DEL", "DUP", "DUP"), 115:118)), pvar)
    write_tsv(tibble(`#IID` = sample_ids), psam)
    scores <- set_names(as.list(sin(1:24)), rev(sample_ids[1:24]))
    write_tsv(bind_cols(tibble(gene_id = gene_ids[[chr]], phenotype_only = 0), as_tibble(scores)), traits)
    .sv_pgen_scan(pgen, pvar, psam, traits, gtf, lof, as.character(chr), max_alt_af = 0.1)
  }) %>% list_rbind()
  genes <- .sv_pgen_gene_results(masks, "BH")
  stopifnot(
    nrow(masks) == 6L, all(masks$status == "ok"), all(masks$n_samples == 24L),
    all(masks$n_lof_carriers == 1L), all(masks$lof_covariate_used),
    identical(masks$n_svs_tested, rep(c(2L, 2L, 4L), 2)),
    nrow(genes) == 2L, all(is.finite(genes$p_value)), all(is.finite(genes$q_value))
  )
  cat("PASS: generated PLINK files, pgenlibr hard calls, SKAT-O, and two-chromosome gathering.\n")
}
run_smoke()
