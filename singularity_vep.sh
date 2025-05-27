#!/bin/bash

INPUT="$1"
if [[ "$INPUT" == *.vcf.gz ]]; then
    BASE_NAME="${INPUT%.vcf.gz}"
elif [[ "$INPUT" == *.vcf ]]; then
    BASE_NAME="${INPUT%.vcf}"
else
    echo "Error: Input file is neither .vcf nor .vcf.gz"
    exit 1
fi

OUTPUT="$(dirname "$INPUT")/$(basename "$BASE_NAME").vep.vcf"

singularity run -H $PWD:/home \
--bind /home/bwubb/resources:/opt/vep/resources \
--bind /home/bwubb/.vep:/opt/vep/.vep \
/appl/containers/vep112.sif vep \
--dir /opt/vep/.vep \
-i $INPUT \
-o $OUTPUT \
--force_overwrite \
--offline \
--cache \
--format vcf \
--vcf --everything --canonical \
--assembly GRCh38 \
--species homo_sapiens \
--fasta /opt/vep/resources/Genomes/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--vcf_info_field ANN \
--plugin NMD \
--plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
--plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
--plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
--custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar_20250106.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
--plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
--plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
