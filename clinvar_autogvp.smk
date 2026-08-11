import os

os.makedirs("data",exist_ok=True)
os.makedirs("results",exist_ok=True)
os.makedirs("work",exist_ok=True)

# Default target: snakemake -s clinvar_autogvp.smk
#   -> work/clinvar_date.txt from data/clinvar.vcf.gz (##fileDate)
# Force refresh: snakemake --forcerun update_clinvar_date
#
# Full build: snakemake -s clinvar_autogvp.smk build_clinvar_autogvp

DATE=""
if os.path.isfile("work/clinvar_date.txt"):
    DATE=open("work/clinvar_date.txt").read().strip()



rule get_clinvar_date:
    input:
        "work/clinvar_date.txt"

rule build_clinvar_autogvp:
    input:
        expand("results/clinvar.autogvp.{date}.vcf.gz",date=DATE),
        expand("results/clinvar.autogvp.{date}.vcf.gz.tbi",date=DATE)

rule download_clinvar_vcf:
    output:
        "data/clinvar.vcf.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
        """

rule update_clinvar_date:
    input:
        "data/clinvar.vcf.gz"
    output:
        "work/clinvar_date.txt"
    shell:
        """
        bcftools view -h {input} | grep "^##fileDate=" | sed 's/##fileDate=//' | sed 's/-//g' > {output}
        """

rule index_clinvar_vcf:
    input:
        "data/clinvar.vcf.gz"
    output:
        "data/clinvar.vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input}
        """

rule download_variant_summary:
    output:
        "data/variant_summary.txt.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
        """

rule download_submission_summary:
    output:
        "data/submission_summary.txt.gz"
    shell:
        """
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz
        """

rule run_clinvar_selection:
    input:
        variant_summary="data/variant_summary.txt.gz",
        submission_summary="data/submission_summary.txt.gz"
    output:
        "results/ClinVar-selected-submissions.tsv"
    params:
        outdir="results"
    shell:
        """
        Rscript ~/software/AutoGVP/scripts/select-clinVar-submissions.R \
        --variant_summary {input.variant_summary} \
        --submission_summary {input.submission_summary} \
        --outdir {params.outdir} \
        --conflict_res latest
        """

rule make_autogvp_tsv:
    input:
        "results/ClinVar-selected-submissions.tsv"
    output:
        "work/autogvp.tsv"
    shell:
        """
        cut -f1,2 {input} > {output}
        """

rule query_clinvar_vcf:
    input:
        vcf="data/clinvar.vcf.gz",
        autogvp="work/autogvp.tsv"
    output:
        "work/clinvar_variants.txt"
    params:
        ids="work/clinvar_ids.txt"
    shell:
        """
        cut -f1 {input.autogvp} | sort -u > {params.ids}
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' -i 'ID=@{params.ids}' {input.vcf} > {output}
        """

rule run_autogvp:
    input:
        autogvp="work/autogvp.tsv",
        clinvar="work/clinvar_variants.txt"
    output:
        "work/clinvar_autogvp.tsv"
    shell:
        """
        declare -A annotation_data
        while IFS=$'\t' read -r clinvar_id annotation; do
            annotation_data["$clinvar_id"]="$annotation"
        done < {input.autogvp}
        while IFS=$'\t' read -r chrom pos clinvar_id ref alt; do
            if [[ -n "${{annotation_data[$clinvar_id]}}" ]]; then
                echo -e "$chrom\t$pos\t$clinvar_id\t$ref\t$alt\t${{annotation_data[$clinvar_id]}}"
            fi
        done < {input.clinvar} > {output}
        """

rule bgzip_autogvp:
    input:
        "work/clinvar_autogvp.tsv"
    output:
        "work/clinvar_autogvp.tsv.gz"
    shell:
        """
        bgzip -c {input} > {output}
        """

rule tabix_autogvp:
    input:
        "work/clinvar_autogvp.tsv.gz"
    output:
        "work/clinvar_autogvp.tsv.gz.tbi"
    shell:
        """
        tabix -s1 -b2 -e2 {input}
        """

rule create_header:
    output:
        "work/header.txt"
    shell:
        """
        echo '##INFO=<ID=AutoGVP,Number=1,Type=String,Description="AutoGVP annotation">' > {output}
        """

rule annotate_clinvar:
    input:
        clinvar_vcf="data/clinvar.vcf.gz",
        clinvar_tbi="data/clinvar.vcf.gz.tbi",
        autogvp="work/clinvar_autogvp.tsv.gz",
        autogvp_tbi="work/clinvar_autogvp.tsv.gz.tbi",
        header="work/header.txt",
        date_file="work/clinvar_date.txt"
    output:
        "results/clinvar.autogvp.{date}.vcf.gz",
        "results/clinvar.autogvp.{date}.vcf.gz.tbi"
    shell:
        """
        bcftools annotate \
        -a {input.autogvp} \
        -c CHROM,POS,ID,REF,ALT,AutoGVP \
        -h {input.header} \
        -Oz -o {output} \
        -W=tbi \
        {input.clinvar_vcf}

        #rsync -vr {output} ~/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz
        #tabix -p vcf ~/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz
        """
