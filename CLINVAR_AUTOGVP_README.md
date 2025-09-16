# ClinVar AutoGVP Pipeline

This Snakemake pipeline updates ClinVar VCF files with AutoGVP annotations.

## Overview

The pipeline performs the following steps:
1. Downloads updated ClinVar VCF, variant_summary, and submission_summary from ClinVar FTP
2. Extracts date from ClinVar VCF header (##fileDate)
3. Downloads AutoGVP annotation file
4. Runs R script to select ClinVar submissions
5. Extracts relevant columns from results
6. Queries ClinVar VCF for variant information
7. Matches ClinVar variants with AutoGVP annotations (integrated into Snakemake)
8. Compresses and indexes AutoGVP results
9. Annotates ClinVar VCF with AutoGVP data

## Prerequisites

- Snakemake
- bcftools
- tabix/bgzip
- R with required packages
- Python 3
- wget

## Configuration

Edit `clinvar_config.yml` to configure:
- Path to concept ID list
- Conflict resolution strategy
- URL for AutoGVP annotation file

## Usage

### Basic run:
```bash
snakemake -s clinvar_autogvp.smk --configfile clinvar_config.yml
```

### Run with cluster:
```bash
snakemake -s clinvar_autogvp.smk --configfile clinvar_config.yml --cluster-config cluster-clinvar.yml --cluster "sbatch --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.cpus}"
```

### Dry run:
```bash
snakemake -s clinvar_autogvp.smk --configfile clinvar_config.yml -n
```

## Output

The final output is a ClinVar VCF file annotated with AutoGVP data:
`results/clinvar/clinvar.autogvp.YYYYMMDD.vcf.gz`

## Files

- `clinvar_autogvp.smk`: Main Snakemake pipeline
- `clinvar_config.yml`: Configuration file
- `cluster-clinvar.yml`: Cluster configuration
- `scripts/select-clinVar-submissions.R`: R script for ClinVar selection

## Notes

- The pipeline extracts the date from the ClinVar VCF header (##fileDate) and uses it for output file naming
- All intermediate files are stored in `work/clinvar/`
- Final results are stored in `results/clinvar/`
- AutoGVP matching is integrated into the Snakemake pipeline using bash associative arrays
- Update the `autogvp_url` in `clinvar_config.yml` with the actual URL for the AutoGVP annotation file
