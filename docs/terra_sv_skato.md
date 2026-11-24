# LoF-adjusted SV SKAT-O on Terra and Researcher Workbench

Use [sv_skato.wdl](../wdl/sv_skato.wdl) with Cromwell. The WDL uses version 1.0 and defaults to `docker.io/epadhi/raresv:workbench-sv-masks-b0cc7fc`. This image includes `R.utils` for gzipped LoF tables. Set `docker_image` if you need a different image.

## Prepare the phenotype BED files

If your expression scores are in a wide TSV, use [the BED splitter](../scripts/split_phenotype_bed.R) to make one `.bed.gz` per chromosome. The expression table needs `gene_id` and one column per participant. The TSS table needs `gene_id`, `chrom`, and `tss_1based`.

```sh
Rscript scripts/split_phenotype_bed.R expression.tsv gene_tss.tsv expression_by_chromosome
```

The splitter matches stable Ensembl IDs, preserves the input gene and sample IDs, and reports expression rows with no TSS match. It writes a one-base BED interval from `tss_1based - 1` to `tss_1based`. Use the TSS table for this preparation step. The Terra workflow reads the phenotype BED directly.

## Chromosome manifest

Prepare a tab-delimited manifest with one row per chromosome. Use these exact headers.

```tsv
chromosome	pgen	pvar	psam	phenotype_bed
1	gs://BUCKET/chr1.pgen	gs://BUCKET/chr1.pvar	gs://BUCKET/chr1.psam	gs://BUCKET/expression_chr1.bed.gz
2	gs://BUCKET/chr2.pgen	gs://BUCKET/chr2.pvar	gs://BUCKET/chr2.psam	gs://BUCKET/expression_chr2.bed.gz
```

- Use one PLINK 2 fileset per chromosome. The `.pvar` must have a `#CHROM` header and preserve the original ALT alleles. Keep `INFO/END` for symbolic deletions and duplications. The pipeline uses the same affected-span rules as the VCF helper. It logs the count of records with unresolved spans.
- The pipeline reads diploid hard calls with `pgenlibr`. Missing hard calls remain missing by default. If every missing hard call represents a confirmed noncarrier, set `missing_pgen_as_reference` to `true` in Workbench or pass `--missing-pgen-as-reference` to the input generator. The scan uses dosage zero for all missing SV hard calls before genotype QC.
- Each phenotype BED starts with `#chr`, `start`, `end`, and `gene_id`, followed by one column per participant and an optional `gene_symbol` column. Use one row per gene and the same chromosome as the PLINK fileset. The sample column names must match `IID` in the PSAM. Plain `.bed` and gzip-compressed `.bed.gz` files are accepted.
- Use the same reference assembly for the phenotype BED and the PVAR. Use that assembly for the GTF if you supply one.
- The shared LoF TSV requires `sample_id`, `gene_id`, and `has_lof_variant`. An omitted gene–participant pair means noncarrier by default. Explicit missing flags are unknown. To model unlisted pairs as unknown, set `lof_missing_as_noncarrier` to `false` in the workflow inputs.
- Optional fixed covariates use a shared TSV with one row per `sample_id`. Set `covariates` and `covariate_columns` together.

If only `./.` and `.|.` encode reference, run [the VCF genotype filter](../scripts/normalize_merged_sv_genotypes.py) before PLINK import. The filter preserves partial calls and single-dot calls. Use `--vcf-half-call missing` so PLINK keeps partial calls missing.

```sh
set -o pipefail
tabix -h input.vcf.gz chr21 \
  | python3 scripts/normalize_merged_sv_genotypes.py \
  | gzip -c > chr21.reference_gt.vcf.gz
plink2 --vcf chr21.reference_gt.vcf.gz --vcf-half-call missing \
  --const-fid 0 --sort-vars --output-chr chrM \
  --make-pgen pvar-cols=xheader,qual,filter,info --out chr21
```

Store the TSV in Cloud Storage. In the workflow form, set `chromosome_manifest` to its `gs://` URI. In Researcher Workbench, enter it in **New job → Prepare inputs**. The manifest task checks the header, checks for duplicate chromosomes, and requires a `gs://` URI in each file column. Cromwell localizes each chromosome's data files as explicit `File` inputs for the scan.

For command-line Cromwell submissions, use [the input generator](../scripts/prepare_sv_skato_inputs.py) to write the inputs JSON. The `chromosome_manifest` value is the URI of the TSV.

```sh
python3 scripts/prepare_sv_skato_inputs.py \
  --manifest gs://BUCKET/chromosomes.tsv \
  --lof-table gs://BUCKET/lof_carriers.tsv \
  --output sv_skato.inputs.json
```

Add `--covariates gs://BUCKET/covariates.tsv --covariate-columns age,PC1,PC2` if needed. Use `--docker-image` to override the default image. Set scan disk and memory to fit the largest chromosome fileset and phenotype BED.

## Model and outputs

By default, the scan extends each BED interval by 1,000,000 bp on both sides and tests SVs whose affected spans overlap that window. Set `window_padding_bp` in Terra or pass `--window-padding-bp` to the manifest converter to change the padding per side. Use `0` for the exact BED interval. The scan limits the left boundary to the start of the chromosome.

To select SVs by CDS overlap, supply `cds_gtf` in Terra or pass `--cds-gtf gs://BUCKET/coding_annotation.gtf.gz` to the manifest converter. The scan then uses the annotated CDS intervals directly. Padding applies only to BED windows.

Each scan matches PSAM IIDs to finite trait scores, then applies the covariate filters. It tests the selected masks with the same quantitative null model per gene. The LoF term enters the model when the flag varies among retained participants. Defaults are alternate-allele frequency ≤1% and missing hard calls ≤5%.

Set `sv_masks` in Workbench or pass `--sv-masks` to the input generator. The default is `DEL,DUP,DEL+DUP`. Separate masks with commas and join types within a mask with `+`. For example, `DEL,DUP,DEL+DUP,INS` adds an insertion-only mask, and `DEL,DUP,INS,DEL+DUP+INS` includes a pooled mask for all three types. The workflow uses the PVAR `POS` as a one-base anchor for each insertion. The PVAR must identify it as an insertion through symbolic `ALT` or `INFO/SVTYPE`. ACAT combines the available p-values from the selected masks, and follow-up uses the same masks.

- `mask_results.tsv` gives raw SKAT-O p-values, sample and SV counts, and mask status. An empty mask has no p-value.
- `gene_results.tsv` gives the equal-weight Cauchy p-value per gene and one FDR correction across all chromosomes. Storey q-values are the default. Set `fdr_method` to `BH` to use Benjamini–Hochberg. If the Storey estimate fails, the task logs the reason and returns missing q-values.
- `significant_genes.tsv` contains genes with `q_value < 0.05`.
- `single_variant_results.tsv` contains one row per gene and SV tested in the follow-up. It gives the SV ID and type, the masks that tested it, the dosage effect and 95% confidence interval, the nominal regression p-value, and the parent gene p-value and q-value.

A gene with one valid mask and no failed masks keeps its mask p-value. A gene with a failed mask has status `incomplete`. It stays in the FDR family at p=1 and has no reported gene p-value or q-value. A gene with no eligible SV test has status `no_test`.

Follow-up runs by default after the gene-level FDR step. Cromwell starts one follow-up task per chromosome with significant genes. The task refits the masks for those genes to recover which SVs SKAT-O used. It fits each gene and SV pair once with the same samples, genotype handling, and covariates, including LoF status when it varies. The variant p-values are nominal because the same data selected the genes. If no genes pass FDR, the combined follow-up table contains only its header. Set `follow_up = false` in Workbench or pass `--no-follow-up` to the input generator to skip this step.

## Verification

```sh
miniwdl check wdl/sv_skato.wdl
python3 tests/test_sv_skato_wdl.py
Rscript tests/test_sv_pgen_workflow.R
Rscript tests/smoke_sv_pgen_workflow.R
```

The [SV SKAT-O image workflow](https://github.com/AoU-Multiomics-Analysis/RareSV/actions/workflows/sv-skato-smoke.yml) builds the pinned micromamba image and runs the WDL and synthetic PLINK checks. Automatic builds run when a script copied into the image changes. Successful builds on `main` publish `ghcr.io/aou-multiomics-analysis/raresv:<commit-sha>`. Pull requests run the checks.

Run the image workflow manually on a branch to test a Dockerfile change. Manual runs on `main` also publish the tested image to GitHub Container Registry. The run summary gives its tag.
