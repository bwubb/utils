#!/bin/bash

INPUT=$1
OUTPUT=`echo $1 | sed -e 's/vcf/vep\.vcf/'`
DIR=`dirname $1`

echo "INPUT: $INPUT"
echo "VEP: $OUTPUT"


#if { conda env list | grep 'vep'; } >/dev/null 2>&1; then conda activate vep; fi

vep -i $INPUT -o $OUTPUT \
--force_overwrite \
--offline \
--cache \
--format vcf \
--vcf \
--everything \
--canonical \
--assembly GRCh38 \
--species homo_sapiens \
--fasta /home/bwubb/resources/Genomes/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--plugin NMD \
--plugin ProteinSeqs,"${DIR}/reference.fa","${DIR}/mutated.fa" \
--plugin Downstream \
--plugin REVEL,/home/bwubb/.vep/revel/revel_grch38.tsv.gz \
--plugin SpliceAI,snv=/home/bwubb/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/home/bwubb/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin gnomADc,/home/bwubb/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
--custom /home/bwubb/.vep/clinvar/vcf_GRCh38/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

bgzip "${OUTPUT}" #&& tabix -fp vcf "${OUTOUT}.gz"

#OUTPUT2=`echo $OUTPUT | sed -e 's/vcf\.gz/report\.csv/'`

#echo "OUTPUT: $OUTPUT2"


#python vep_vcf_parser.py -i $OUTPUT -o $OUTPUT2 everything
