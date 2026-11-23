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
