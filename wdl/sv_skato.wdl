version 1.0

struct ChromosomeInputs {
  String chromosome
  File pgen
  File pvar
  File psam
  File molecular_traits
}

workflow sv_skato {
  input {
    Array[ChromosomeInputs]+ chromosome_manifest
    File cds_gtf
    File lof_table
    File? covariates
    String covariate_columns = ""
    Boolean lof_missing_as_noncarrier = true
    Float max_alt_af = 0.01
    Float max_missing = 0.05
    String fdr_method = "storey"
    String docker_image
    String scan_memory = "12 GB"
    Int scan_disk_gb = 100
    String gather_memory = "4 GB"
    Int gather_disk_gb = 20
  }

  scatter (entry in chromosome_manifest) {
    call test_chromosome {
      input:
        files = entry,
        cds_gtf = cds_gtf,
        lof_table = lof_table,
        covariates = covariates,
        covariate_columns = covariate_columns,
        lof_missing_as_noncarrier = lof_missing_as_noncarrier,
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

task test_chromosome {
  input {
    ChromosomeInputs files
    File cds_gtf
    File lof_table
    File? covariates
    String covariate_columns
    Boolean lof_missing_as_noncarrier
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
    covariate_file='~{if defined(covariates) then sub(select_first([covariates]), "'", "'\"'\"'") else ""}'
    covariate_names='~{sub(covariate_columns, "'", "'\"'\"'")}'
    extra_args=()
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
      --traits '~{sub(files.molecular_traits, "'", "'\"'\"'")}' \
      --gtf '~{sub(cds_gtf, "'", "'\"'\"'")}' \
      --lof '~{sub(lof_table, "'", "'\"'\"'")}' \
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
