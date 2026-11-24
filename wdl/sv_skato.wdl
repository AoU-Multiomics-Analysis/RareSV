version 1.0

struct ChromosomeInputs {
  String chromosome
  File pgen
  File pvar
  File psam
  File phenotype_bed
}

workflow sv_skato {
  input {
    File chromosome_manifest
    File? cds_gtf
    File lof_table
    File? covariates
    String covariate_columns = ""
    Boolean lof_missing_as_noncarrier = true
    Int window_padding_bp = 1000000
    Float max_alt_af = 0.01
    Float max_missing = 0.05
    String fdr_method = "storey"
    String docker_image = "docker.io/epadhi/raresv:olink-base-1445dab"
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
    ChromosomeInputs entry = object {
      chromosome: row[0],
      pgen: pgen_file,
      pvar: pvar_file,
      psam: psam_file,
      phenotype_bed: phenotype_bed_file
    }

    call test_chromosome {
      input:
        files = entry,
        cds_gtf = cds_gtf,
        lof_table = lof_table,
        covariates = covariates,
        covariate_columns = covariate_columns,
        lof_missing_as_noncarrier = lof_missing_as_noncarrier,
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
      masks = test_chromosome.mask_results,
      fdr_method = fdr_method,
      docker_image = docker_image,
      task_memory = gather_memory,
      disk_gb = gather_disk_gb
  }

  output {
    Array[File] chromosome_mask_results = test_chromosome.mask_results
    File mask_results = gather_results.mask_results
    File gene_results = gather_results.gene_results
    File significant_genes = gather_results.significant_genes
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
    ChromosomeInputs files
    File? cds_gtf
    File lof_table
    File? covariates
    String covariate_columns
    Boolean lof_missing_as_noncarrier
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
    chromosome='~{sub(files.chromosome, "'", "'\"'\"'")}'
    cds_annotation='~{if defined(cds_gtf) then sub(select_first([cds_gtf]), "'", "'\"'\"'") else ""}'
    covariate_file='~{if defined(covariates) then sub(select_first([covariates]), "'", "'\"'\"'") else ""}'
    covariate_names='~{sub(covariate_columns, "'", "'\"'\"'")}'
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

    printf 'Starting LoF-adjusted SKAT-O for chromosome %s\n' "$chromosome" >&2
    Rscript /opt/sv-skato/sv_pgen_workflow.R scan \
      --chromosome "$chromosome" \
      --pgen '~{sub(files.pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(files.pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(files.psam, "'", "'\"'\"'")}' \
      --phenotype-bed '~{sub(files.phenotype_bed, "'", "'\"'\"'")}' \
      --lof '~{sub(lof_table, "'", "'\"'\"'")}' \
      --window-padding-bp '~{window_padding_bp}' \
      --max-alt-af '~{max_alt_af}' --max-missing '~{max_missing}' \
      "${extra_args[@]}"
    printf 'Finished chromosome %s\n' "$chromosome" >&2
  >>>

  output {
    File mask_results = "mask_results.tsv"
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
    String fdr_method
    String docker_image
    String task_memory
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    # The task writes this newline-delimited list from localized File inputs during command rendering.
    mask_list='~{sub(write_lines(masks), "'", "'\"'\"'")}'
    printf 'Combining chromosome mask results and applying 5%% FDR across genes\n' >&2
    Rscript /opt/sv-skato/sv_pgen_workflow.R gather \
      --mask-list "$mask_list" --fdr-method '~{sub(fdr_method, "'", "'\"'\"'")}'
    printf 'Finished writing combined gene results\n' >&2
  >>>

  output {
    File mask_results = "mask_results.tsv"
    File gene_results = "gene_results.tsv"
    File significant_genes = "significant_genes.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: task_memory
    disks: "local-disk ~{disk_gb} SSD"
  }
}
