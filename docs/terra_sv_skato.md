# LoF-adjusted SV SKAT-O on Terra

Use [sv_skato.wdl](../wdl/sv_skato.wdl) with the image from [the Dockerfile](../docker/sv-skato/Dockerfile). The WDL uses version 1.0. Supply an image tag or digest that Terra can pull as `docker_image`.

## Chromosome manifest

Prepare a tab-delimited manifest with one row per chromosome. Use these exact headers.

```tsv
chromosome	pgen	pvar	psam	molecular_traits
1	gs://BUCKET/chr1.pgen	gs://BUCKET/chr1.pvar	gs://BUCKET/chr1.psam	gs://BUCKET/chr1.expression.tsv
2	gs://BUCKET/chr2.pgen	gs://BUCKET/chr2.pvar	gs://BUCKET/chr2.psam	gs://BUCKET/chr2.expression.tsv
```

- Use one PLINK 2 fileset per chromosome. The `.pvar` must have a `#CHROM` header and preserve the original ALT alleles. Keep `INFO/END` for symbolic deletions and duplications. The pipeline uses the same affected-span rules as the VCF helper. It logs the count of records with unresolved spans.
- The pipeline reads diploid hard calls with `pgenlibr`. A missing PGEN hard call remains missing. Encode calls that mean confirmed reference as `0/0` before writing the PGEN.
- Each molecular trait TSV has one row per gene, a `gene_id` column, and one column per participant. The sample column names must match `IID` in the PSAM. An optional `gene_symbol` column supplies the display label. The GTF and PVAR must use the same reference assembly.
- The shared LoF TSV requires `sample_id`, `gene_id`, and `has_lof_variant`. An omitted gene–participant pair means noncarrier by default. Explicit missing flags are unknown. To model unlisted pairs as unknown, set `lof_missing_as_noncarrier` to `false` in the workflow inputs.
- Optional fixed covariates use a shared TSV with one row per `sample_id`. Set `covariates` and `covariate_columns` together.

Convert the TSV to workflow inputs with [the manifest converter](../scripts/prepare_sv_skato_inputs.py). The converter creates a typed array for the WDL. Terra localizes the files in each row as explicit `File` inputs.

```sh
python3 scripts/prepare_sv_skato_inputs.py \
  --manifest chromosomes.tsv \
  --cds-gtf gs://BUCKET/coding_annotation.gtf.gz \
  --lof-table gs://BUCKET/lof_carriers.tsv \
  --docker-image REGISTRY/sv-skato:TAG \
  --output sv_skato.inputs.json
```

Add `--covariates gs://BUCKET/covariates.tsv --covariate-columns age,PC1,PC2` if needed. Use the resulting values when you configure the workflow in Terra. Set scan disk and memory to fit the largest chromosome fileset and trait matrix.

## Model and outputs

Each scan matches PSAM IIDs to finite trait scores, then applies the covariate filters. It tests CDS-overlapping `DEL`, `DUP`, and `DEL+DUP` masks with the same quantitative null model per gene. The LoF term enters the model when the flag varies among retained participants. Defaults are alternate-allele frequency ≤1% and missing hard calls ≤5%.

- `mask_results.tsv` gives raw SKAT-O p-values, sample and SV counts, and mask status. An empty mask has no p-value.
- `gene_results.tsv` gives the equal-weight Cauchy p-value per gene and one FDR correction across all chromosomes. Storey q-values are the default. Set `fdr_method` to `BH` to use Benjamini–Hochberg. If the Storey estimate fails, the task logs the reason and returns missing q-values.
- `significant_genes.tsv` contains genes with `q_value < 0.05`.

A gene with one valid mask and no failed masks keeps its mask p-value. A gene with a failed mask has status `incomplete`. It stays in the FDR family at p=1 and has no reported gene p-value or q-value. A gene with no eligible SV test has status `no_test`.

## Verification

```sh
miniwdl check wdl/sv_skato.wdl
python3 tests/test_sv_skato_wdl.py
Rscript tests/test_sv_pgen_workflow.R
```

The [SV SKAT-O image workflow](https://github.com/AoU-Multiomics-Analysis/RareSV/actions/workflows/sv-skato-smoke.yml) builds the pinned micromamba image and runs the WDL and synthetic PLINK checks. Automatic builds run when a script copied into the image changes. Successful builds on `main` publish `ghcr.io/aou-multiomics-analysis/raresv:<commit-sha>`. Pull requests run the checks.

For the first image, run the workflow manually on `main`. The run summary gives the image tag to use for `docker_image`. GitHub Container Registry packages are private on first publication. Use a registry reference that your Terra workspace has permission to pull.
