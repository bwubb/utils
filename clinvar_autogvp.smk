import os
import yaml
from datetime import datetime

### FUNCTIONS ###

def get_clinvar_date():
    """Get current date for ClinVar file naming"""
    return datetime.now().strftime('%Y%m%d')

### ### CONFIGURATION ### ###

# Load configuration
with open(config['clinvar']['config_file'], 'r') as f:
    CLINVAR_CONFIG = yaml.safe_load(f)

# Set up directories
os.makedirs("logs/cluster/clinvar", exist_ok=True)
os.makedirs("data/clinvar", exist_ok=True)
os.makedirs("results/clinvar", exist_ok=True)
os.makedirs("work/clinvar", exist_ok=True)

### ### RULES ### ###

localrules: clinvar_autogvp_all

rule clinvar_autogvp_all:
    input:
        "results/clinvar/clinvar.autogvp.{date}.vcf.gz"

rule download_clinvar_vcf:
    output:
        "data/clinvar/clinvar.vcf.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
        """

rule download_variant_summary:
    output:
        "data/clinvar/variant_summary.txt.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
        """

rule download_submission_summary:
    output:
        "data/clinvar/submission_summary.txt.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
        """

rule run_clinvar_selection:
    input:
        variant_summary="data/clinvar/variant_summary.txt.gz",
        submission_summary="data/clinvar/submission_summary.txt.gz",
        concept_ids=CLINVAR_CONFIG['concept_id_list']
    output:
        "results/clinvar/selected_clinvar_submissions.txt"
    params:
        outdir="results/clinvar",
        conflict_res=CLINVAR_CONFIG.get('conflict_resolution', 'latest')
    shell:
        """
        Rscript scripts/select-clinVar-submissions.R \
        --variant_summary {input.variant_summary} \
        --submission_summary {input.submission_summary} \
        --outdir {params.outdir} \
        --conceptID_list {input.concept_ids} \
        --conflict_res {params.conflict_res}
        """

rule extract_columns:
    input:
        "results/clinvar/selected_clinvar_submissions.txt"
    output:
        "work/clinvar/selected_columns.txt"
    shell:
        """
        cut -f1,2 {input} > {output}
        """

rule query_clinvar_vcf:
    input:
        vcf="data/clinvar/clinvar.vcf.gz",
        columns="work/clinvar/selected_columns.txt"
    output:
        "work/clinvar/clinvar_variants.txt"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {input.vcf} > {output}
        """

rule extract_clinvar_date:
    input:
        "data/clinvar/clinvar.vcf.gz"
    output:
        "work/clinvar/clinvar_date.txt"
    shell:
        """
        bcftools view -h {input} | grep "^##fileDate=" | sed 's/##fileDate=//' | sed 's/-//g' > {output}
        """

rule download_autogvp:
    output:
        "work/clinvar/autogvp.tsv"
    params:
        url=CLINVAR_CONFIG.get('autogvp_url', 'https://example.com/autogvp.tsv')
    shell:
        """
        wget -O {output} {params.url}
        """

rule run_autogvp:
    input:
        autogvp="work/clinvar/autogvp.tsv",
        clinvar="work/clinvar/clinvar_variants.txt"
    output:
        "work/clinvar/clinvar_autogvp.tsv"
    shell:
        """
        # Load AutoGVP annotations into dictionary
        declare -A annotation_data
        while IFS=$'\t' read -r clinvar_id annotation; do
            annotation_data["$clinvar_id"]="$annotation"
        done < {input.autogvp}
        
        # Process ClinVar variants and match with AutoGVP
        while IFS=$'\t' read -r chrom pos clinvar_id ref alt; do
            if [[ -n "${{annotation_data[$clinvar_id]}}" ]]; then
                echo -e "$chrom\t$pos\t$clinvar_id\t$ref\t$alt\t${{annotation_data[$clinvar_id]}}"
            fi
        done < {input.clinvar} > {output}
        """

rule bgzip_autogvp:
    input:
        "work/clinvar/clinvar_autogvp.tsv"
    output:
        "work/clinvar/clinvar_autogvp.tsv.gz"
    shell:
        """
        bgzip -c {input} > {output}
        """

rule tabix_autogvp:
    input:
        "work/clinvar/clinvar_autogvp.tsv.gz"
    output:
        "work/clinvar/clinvar_autogvp.tsv.gz.tbi"
    shell:
        """
        tabix -s1 -b2 -e2 {input}
        """

rule create_header:
    output:
        "work/clinvar/header.txt"
    shell:
        """
        echo '##INFO=<ID=AutoGVP,Number=1,Type=String,Description="AutoGVP annotation">' > {output}
        """

rule annotate_clinvar:
    input:
        vcf="data/clinvar/clinvar.vcf.gz",
        autogvp="work/clinvar/clinvar_autogvp.tsv.gz",
        header="work/clinvar/header.txt",
        date_file="work/clinvar/clinvar_date.txt"
    output:
        "results/clinvar/clinvar.autogvp.{date}.vcf.gz"
    params:
        date=lambda wildcards, input: open(input.date_file).read().strip()
    shell:
        """
        bcftools annotate \
        -a {input.autogvp} \
        -c CHROM,POS,ID,REF,ALT,AutoGVP \
        -h {input.header} \
        -Oz -o {output} \
        {input.vcf}
        """

rule index_final_vcf:
    input:
        "results/clinvar/clinvar.autogvp.{date}.vcf.gz"
    output:
        "results/clinvar/clinvar.autogvp.{date}.vcf.gz.tbi"
    shell:
        """
        bcftools index {input}
        """
