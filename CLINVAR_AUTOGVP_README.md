# ClinVar AutoGVP Pipeline

This Snakemake pipeline builds an AutoGVP-annotated ClinVar VCF for downstream use (for example VEP and variant parsing that read the `AutoGVP` INFO field). AutoGVP labels here are not downloaded from a separate resource: they come from the ClinVar submission selection step, which applies your gene/concept list and conflict-resolution rules to ClinVar’s tab-delimited releases.

## Overview

The pipeline performs the following steps:

1. Downloads the current ClinVar GRCh38 VCF, `variant_summary.txt.gz`, and `submission_summary.txt.gz` from ClinVar FTP into `data/clinvar/`.
2. Runs `scripts/select-clinVar-submissions.R` on those tables with the configured concept ID list and conflict resolution. The R script writes **`results/ClinVar-selected-submissions.tsv`** (under `--outdir results`). That file is the authoritative selection table for this workflow.
3. Builds **`work/clinvar/autogvp.tsv`** by taking the **first two tab-separated columns** of the selection TSV. Column 1 is the ClinVar variation ID; column 2 is the AutoGVP annotation string used later in the VCF INFO field. This replaces the old `extract_columns` / `selected_columns.txt` path and is the AutoGVP lookup table for the rest of the pipeline.
4. Extracts **`##fileDate`** from the ClinVar VCF header and normalizes it (no hyphens) for naming the final output (`YYYYMMDD`).
5. Queries the ClinVar VCF for rows whose **`ID`** appears in column 1 of `autogvp.tsv`, writing **`work/clinvar/clinvar_variants.txt`** (CHROM, POS, ID, REF, ALT).
6. Joins variant coordinates with AutoGVP strings in **`work/clinvar/clinvar_autogvp.tsv`** (bash associative array keyed by ClinVar ID).
7. Compresses and indexes that annotation file (`bgzip` + `tabix`) for `bcftools annotate`.
8. Adds an `AutoGVP` INFO definition and annotates the full ClinVar VCF, producing the dated result under **`results/clinvar/`**.

## Prerequisites

- Snakemake
- bcftools
- tabix / bgzip
- R with packages required by `select-clinVar-submissions.R`
- wget
- bash (associative arrays for the join in `run_autogvp`)

## Configuration

Edit `clinvar_config.yml` and pass it with `--configfile`. Snakemake loads it as `config`; the workflow reads the `clinvar` section only.

- **`concept_id_list`**: Path to the concept ID list passed to the R script (`--conceptID_list`).
- **`conflict_resolution`**: How to resolve conflicting submissions when multiple apply (default in the workflow: `latest`). Passed to the R script as `--conflict_res`.

The `ftp` and `output` blocks in the config file document naming conventions; ClinVar download URLs are defined in `clinvar_autogvp.smk`.

## Usage

### Basic run

```bash
snakemake -s clinvar_autogvp.smk --configfile clinvar_config.yml
```

### Run on a cluster

```bash
snakemake -s clinvar_autogvp.smk --configfile clinvar_config.yml \
  --cluster-config cluster-clinvar.yml \
  --cluster "sbatch --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.cpus}"
```

### Dry run

```bash
snakemake -s clinvar_autogvp.smk --configfile clinvar_config.yml -n
```

## Output

**Primary output:** ClinVar VCF with `AutoGVP` in INFO:

`results/clinvar/clinvar.autogvp.YYYYMMDD.vcf.gz`

(and `.tbi` from `index_final_vcf`)

**Key intermediates:**

| Path | Role |
|------|------|
| `results/ClinVar-selected-submissions.tsv` | R selection output (must match R script name and `--outdir results`) |
| `work/clinvar/autogvp.tsv` | Two-column ID + AutoGVP annotation lookup |
| `work/clinvar/clinvar_variants.txt` | Subset of ClinVar VCF rows for selected IDs |
| `work/clinvar/clinvar_autogvp.tsv.gz` | Annotation file for `bcftools annotate` |

## Files

- `clinvar_autogvp.smk` — main workflow
- `clinvar_config.yml` — concept list and conflict resolution
- `cluster-clinvar.yml` — optional cluster resources
- `scripts/select-clinVar-submissions.R` — ClinVar submission selection (must write `ClinVar-selected-submissions.tsv` into the configured outdir)

## Notes

- The final VCF date suffix comes from the ClinVar release **`##fileDate`** in `data/clinvar/clinvar.vcf.gz`, not from the run date on the cluster.
- Cached downloads and working files live under `data/clinvar/` and `work/clinvar/`; selection TSV stays in `results/` because that is where the R script writes it.
- If the R script’s output name or outdir changes, update **`run_clinvar_selection`** and **`make_autogvp_tsv`** inputs/outputs in `clinvar_autogvp.smk` to match—Snakemake will not guess the filename.
- Only variants present in both the ClinVar VCF and `autogvp.tsv` receive an `AutoGVP` INFO value; others are left unchanged by the annotate step.
