#needs plan_mnp_gt.py
#need find_mnp_variants.py
#needs test_mnp_sample_info.py


rule all:
    input:
        "data/mnp/mnp_read_check.assignments.tsv",
        "data/mnp/mnp_read_check.samples.txt",
        #expand("data/mnp/chr{CHR}.mnp_gt.plan.txt",CHR=range(1,23))



rule find_mnp_variants:
    input:
        csv="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.report.csv"
    output:
        csv="data/mnp/chr{CHR}.annotation.mnp.csv",
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    shell:
        "python find_mnp_variants.py -i {input.csv} -o {output.csv} -m no_sample -r {output.regions} -I {output.id_file}"

rule bcftools_mnp_filter:
    input:
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.bcf",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    output:
        vcf="data/mnp/chr{CHR}.mnp.vcf"
    shell:
        """
        bcftools view -R {input.regions} -i 'ID=@{input.id_file}' -Ov -o {output.vcf} {input.vcf}
        """

rule bcftools_mnp_sample_info:
    input:
        vcf="data/mnp/chr{CHR}.mnp.vcf"
    output:
        "data/mnp/chr{CHR}.mnp_sample_info.txt"
    shell:
        """
        bcftools query -f '[%ID %SAMPLE %GT %AD %DP %VAF\n]' {input} |
        awk -F' ' '$3=="0/1" || $3=="1/1"' > {output}
        """

rule test_mnp_sample_info:
    input:
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.txt",
    output:
        pass_out="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",
        fail_out="data/mnp/chr{CHR}.mnp_sample_info.FAIL.txt"
    params:
        case_list="case.list"
    shell:
        "python test_mnp_sample_info.py -p {input.pairs} -i {input.sample_info} -o {output.pass_out} -f {output.fail_out} --case-list {params.case_list}"


rule choose_mnp_read_check_samples:
    input:
        pass_files=expand("data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",CHR=range(1,23)),
    output:
        assignments="data/mnp/mnp_read_check.assignments.tsv",
        samples="data/mnp/mnp_read_check.samples.txt",
    shell:
        "python choose_mnp_read_check_samples.py {input.pass_files} -o {output.assignments} -s {output.samples}"

#rule plan_mnp_gt:
#    input:
#        sample_info="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt"
#    output:
#        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt"
#    params:
#        ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
#    shell:
#        "python plan_mnp_gt.py -i {input.sample_info} -o {output.plan} -r {params.ref}"

#rule manage_mnp_gt:
#    input:
#        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt",
#        vcf="data/mnp/chr{CHR}.mnp.vcf"
#    output:
#        vcf="data/mnp/chr{CHR}.mnp_gt.vcf"
#    shell:
#        "python manage_mnp_gt.py -p {input.plan} -i {input.vcf} -o {output.vcf}"
#
#rule annotate_mnp_gt:
#    input:
#        vcf="data/mnp/chr{CHR}.mnp_gt.vcf"
#    output:
#        vcf="data/mnp/chr{CHR}.mnp_gt.vep.vcf"
#    shell:
#        """
#        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
#        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
#        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
#        --bind /home/bwubb/resources:/opt/vep/resources \
#        --bind /home/bwubb/.vep:/opt/vep/.vep \
#        /appl/containers/vep112.sif vep \
#        --dir /opt/vep/.vep \
#        -i $PWD/{input} \
#        -o $PWD/{output} \
#        --force_overwrite \
#        --offline \
#        --cache \
#        --format vcf \
#        --vcf --everything --canonical \
#        --assembly GRCh38 \
#        --species homo_sapiens \
#        --fasta /opt/vep/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta \
#        --vcf_info_field ANN \
#        --plugin NMD \
#        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
#        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
#        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
#        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
#        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
#        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
#        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
#        """
#
#rule parse_mnp_gt_vep:
#    input:
#        vcf="data/mnp/chr{CHR}.mnp_gt.vep.vcf"
#    output:
#        csv="data/mnp/chr{CHR}.mnp_gt.vep.csv"
#    shell:
#        "python vep_vcf_parser2.py -i {input.vcf} -o {output.csv} -m cohort"
  