version 1.0

workflow sv_skato {
  input {
    File chromosome_manifest
    File? cds_gtf
    File lof_table
    File? covariates
    String covariate_columns = ""
    Boolean lof_missing_as_noncarrier = true
    Boolean missing_pgen_as_reference = false
    String sv_masks = "DEL,DUP,DEL+DUP"
    Int window_padding_bp = 1000000
    Float max_alt_af = 0.01
    Float max_missing = 0.05
    String fdr_method = "storey"
    Boolean follow_up = true
    String docker_image = "docker.io/epadhi/raresv:workbench-sv-masks-b0cc7fc"
    String scan_memory = "12 GB"
    Int scan_disk_gb = 100
    String gather_memory = "4 GB"
    Int gather_disk_gb = 20
  }

  call parse_manifest {
    input:
      manifest = chromosome_manifest,
      docker_image = docker_image
  }

  scatter (row in parse_manifest.rows) {
    File pgen_file = row[1]
    File pvar_file = row[2]
    File psam_file = row[3]
    File phenotype_bed_file = row[4]
    call test_chromosome {
      input:
        chromosome = row[0],
        pgen = pgen_file,
        pvar = pvar_file,
        psam = psam_file,
        phenotype_bed = phenotype_bed_file,
        cds_gtf = cds_gtf,
        lof_table = lof_table,
        covariates = covariates,
        covariate_columns = covariate_columns,
        lof_missing_as_noncarrier = lof_missing_as_noncarrier,
        missing_pgen_as_reference = missing_pgen_as_reference,
        sv_masks = sv_masks,
        window_padding_bp = window_padding_bp,
        max_alt_af = max_alt_af,
        max_missing = max_missing,
        docker_image = docker_image,
        task_memory = scan_memory,
        disk_gb = scan_disk_gb
    }
  }

  call gather_results {
    input:
      masks = test_chromosome.results,
      manifest = parse_manifest.parsed_manifest,
      sv_masks = sv_masks,
      fdr_method = fdr_method,
      docker_image = docker_image,
      task_memory = gather_memory,
      disk_gb = gather_disk_gb
  }

  if (follow_up) {
    scatter (follow_up_row in gather_results.follow_up_rows) {
      File follow_up_pgen = follow_up_row[1]
      File follow_up_pvar = follow_up_row[2]
      File follow_up_psam = follow_up_row[3]
      File follow_up_bed = follow_up_row[4]
      call test_chromosome as follow_up_chromosome {
        input:
          chromosome = follow_up_row[0],
          pgen = follow_up_pgen,
          pvar = follow_up_pvar,
          psam = follow_up_psam,
          phenotype_bed = follow_up_bed,
          cds_gtf = cds_gtf,
          lof_table = lof_table,
          covariates = covariates,
          covariate_columns = covariate_columns,
          lof_missing_as_noncarrier = lof_missing_as_noncarrier,
          missing_pgen_as_reference = missing_pgen_as_reference,
          sv_masks = sv_masks,
          window_padding_bp = window_padding_bp,
          max_alt_af = max_alt_af,
          max_missing = max_missing,
          significant_genes = gather_results.significant_genes,
          docker_image = docker_image,
          task_memory = scan_memory,
          disk_gb = scan_disk_gb
      }
    }
    call gather_follow_up {
      input:
        files = follow_up_chromosome.results,
        docker_image = docker_image,
        task_memory = gather_memory,
        disk_gb = gather_disk_gb
    }
  }

  output {
    Array[File] chromosome_mask_results = test_chromosome.results
    File mask_results = gather_results.mask_results
    File gene_results = gather_results.gene_results
    File significant_genes = gather_results.significant_genes
    File? single_variant_results = gather_follow_up.results
  }
}

task parse_manifest {
  input {
    File manifest
    String docker_image
  }

  command <<<
    set -euo pipefail
    manifest_file='~{sub(manifest, "'", "'\"'\"'")}'
    printf 'Reading chromosome manifest\n' >&2
    if [[ "$manifest_file" == *://* || ! -f "$manifest_file" || ! -r "$manifest_file" ]]; then
      printf 'Localization error for chromosome manifest: file is not readable: %s\n' "$manifest_file" >&2
      exit 1
    fi

    LC_ALL=C awk '
      BEGIN {
        FS = OFS = "\t"
        split("chromosome pgen pvar psam phenotype_bed", columns, " ")
      }
      function fail(detail) {
        print "Invalid chromosome manifest: " detail > "/dev/stderr"
        failed = 1
        exit 1
      }
      { sub(/\r$/, "") }
      NR == 1 {
        sub(/^\357\273\277/, "", $1)
        if (NF != 5) fail("expected chromosome, pgen, pvar, psam, and phenotype_bed columns")
        for (i = 1; i <= NF; i++) index_of[$i] = i
        for (i = 1; i <= 5; i++) {
          if (!(columns[i] in index_of)) fail("missing column " columns[i])
        }
        next
      }
      /^[[:space:]]*$/ { next }
      {
        if (NF != 5 || /\r/) fail("expected five tab-delimited fields at line " NR)
        for (i = 1; i <= 5; i++) {
          value[i] = $(index_of[columns[i]])
          sub(/^[[:space:]]+/, "", value[i])
          sub(/[[:space:]]+$/, "", value[i])
        }
        key = toupper(value[1])
        sub(/^CHR/, "", key)
        if (key == "23") key = "X"
        if (key == "24") key = "Y"
        if (key == "26" || key == "M") key = "MT"
        if (key == "" || seen[key]++) fail("blank or repeated chromosome at line " NR)
        for (i = 2; i <= 5; i++) {
          if (value[i] !~ /^gs:\/\/[^\/[:space:]]+\/.+$/) {
            fail("expected a gs:// path for " columns[i] " at line " NR)
          }
        }
        print value[1], value[2], value[3], value[4], value[5]
        rows++
      }
      END { if (!failed && !rows) fail("expected at least one chromosome row") }
    ' < "$manifest_file" > manifest_rows.tsv
    printf 'Validated %s chromosome rows\n' "$(wc -l < manifest_rows.tsv)" >&2
  >>>

  output {
    Array[Array[String]] rows = read_tsv("manifest_rows.tsv")
    File parsed_manifest = "manifest_rows.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 SSD"
  }
}

task test_chromosome {
  input {
    String chromosome
    File pgen
    File pvar
    File psam
    File phenotype_bed
    File? cds_gtf
    File lof_table
    File? covariates
    File? significant_genes
    String covariate_columns
    Boolean lof_missing_as_noncarrier
    Boolean missing_pgen_as_reference
    String sv_masks
    Int window_padding_bp
    Float max_alt_af
    Float max_missing
    String docker_image
    String task_memory
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
    chromosome='~{sub(chromosome, "'", "'\"'\"'")}'
    cds_annotation='~{if defined(cds_gtf) then sub(select_first([cds_gtf]), "'", "'\"'\"'") else ""}'
    covariate_file='~{if defined(covariates) then sub(select_first([covariates]), "'", "'\"'\"'") else ""}'
    covariate_names='~{sub(covariate_columns, "'", "'\"'\"'")}'
    follow_up_file='~{if defined(significant_genes) then sub(select_first([significant_genes]), "'", "'\"'\"'") else ""}'
    extra_args=()
    if [[ -n "$cds_annotation" ]]; then
      extra_args+=(--gtf "$cds_annotation")
    fi
    if [[ -n "$covariate_file" ]]; then
      extra_args+=(--covariates "$covariate_file")
    fi
    if [[ -n "$covariate_names" ]]; then
      extra_args+=(--covariate-cols "$covariate_names")
    fi
    if [[ '~{lof_missing_as_noncarrier}' == false ]]; then
      extra_args+=(--unknown-unlisted-lof)
    fi
    if [[ '~{missing_pgen_as_reference}' == true ]]; then
      extra_args+=(--missing-pgen-as-reference)
    fi

    if [[ -n "$follow_up_file" ]]; then
      extra_args+=(--significant-genes "$follow_up_file")
      printf 'Starting single-variant follow-up for chromosome %s\n' "$chromosome" >&2
    else
      printf 'Starting LoF-adjusted SKAT-O for chromosome %s\n' "$chromosome" >&2
    fi
    Rscript /opt/sv-skato/sv_pgen_workflow.R scan \
      --chromosome "$chromosome" \
      --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(psam, "'", "'\"'\"'")}' \
      --phenotype-bed '~{sub(phenotype_bed, "'", "'\"'\"'")}' \
      --lof '~{sub(lof_table, "'", "'\"'\"'")}' \
      --sv-masks '~{sub(sv_masks, "'", "'\"'\"'")}' \
      --window-padding-bp '~{window_padding_bp}' \
      --max-alt-af '~{max_alt_af}' --max-missing '~{max_missing}' \
      "${extra_args[@]}"
    printf 'Finished chromosome %s\n' "$chromosome" >&2
  >>>

  output {
    File results = if defined(significant_genes) then "single_variant_results.tsv" else "mask_results.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: task_memory
    disks: "local-disk ~{disk_gb} SSD"
  }
}

task gather_results {
  input {
    Array[File] masks
    File manifest
    String fdr_method
    String sv_masks
    String docker_image
    String task_memory
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    # Pass paths as literal stdin. Replace the outer shell so it cannot run anything after the delimiter.
    exec /bin/bash -euo pipefail -c '
      mask_list="mask_paths.txt"
      cat > "$mask_list"
      if [[ $(wc -l < "$mask_list") -ne "$2" ]]; then
        printf "Expected one localized mask path per input file\n" >&2
        exit 1
      fi
      printf "Combining chromosome mask results and applying 5%% FDR across genes\n" >&2
      Rscript /opt/sv-skato/sv_pgen_workflow.R gather --mask-list "$mask_list" --manifest "$3" --sv-masks "$4" --fdr-method "$1"
      printf "Finished writing combined gene results\n" >&2
    ' gather-results '~{sub(fdr_method, "'", "'\"'\"'")}' '~{length(masks)}' '~{sub(manifest, "'", "'\"'\"'")}' '~{sub(sv_masks, "'", "'\"'\"'")}' <<'RARESV_MASK_PATHS'
~{sep="\n" masks}
RARESV_MASK_PATHS
  >>>

  output {
    File mask_results = "mask_results.tsv"
    File gene_results = "gene_results.tsv"
    File significant_genes = "significant_genes.tsv"
    Array[Array[String]] follow_up_rows = if size("follow_up_manifest.tsv", "B") > 0.0 then read_tsv("follow_up_manifest.tsv") else []
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: task_memory
    disks: "local-disk ~{disk_gb} SSD"
  }
}

task gather_follow_up {
  input {
    Array[File] files
    String docker_image
    String task_memory
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    # Read the rendered paths as literal input and keep the list in the task directory.
    exec /bin/bash -euo pipefail -c '
      file_list="follow_up_paths.txt"
      if [[ "$1" -eq 0 ]]; then
        : > "$file_list"
      else
        cat > "$file_list"
        if [[ $(wc -l < "$file_list") -ne "$1" ]]; then
          printf "Expected one localized follow-up path per input file\n" >&2
          exit 1
        fi
      fi
      printf "Combining single-variant follow-up results\n" >&2
      Rscript /opt/sv-skato/sv_pgen_workflow.R gather-follow-up --follow-up-list "$file_list"
    ' gather-follow-up '~{length(files)}' <<'RARESV_FOLLOWUP_PATHS'
~{sep="\n" files}
RARESV_FOLLOWUP_PATHS
  >>>

  output {
    File results = "single_variant_results.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "~{task_memory}"
    disks: "local-disk ~{disk_gb} SSD"
  }
}
