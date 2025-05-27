include: "common_snps.smk"
include: "metrics.smk"

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()

with open(config.get('project',{}).get('bam_table','bams.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

id=config['resources']['targets_key']

rule collect_ibd:
    input:
        "data/work/IBD/ibd-related.txt"

rule gatk_alleles:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample],
        snp=f"data/work/common_snps/{id}/gnomad.exomes.v4.sites.{id}.common_biallelic_snps.simplified.vcf.gz"
    output:
        "data/work/{sample}/gatk/common_snps.g.vcf.gz"
    params:
        ref=config['referemce']['fasta']
    shell:
        """
        gatk --java-options '-Xmx10240m' HaplotypeCaller \
        -R {params.ref} \
        -I {input.bam} \
        -O {output} \
        -L {input.snp} \
        --alleles {input.snp} \
        --genotype-filtered-alleles true \
        --emit-ref-confidence GVCF \
        --output-mode EMIT_ALL_CONFIDENT_SITES \
        --interval-padding 0
        """

rule make_sample_map:
    input:
        expand("data/work/{sample}/gatk/common_snps.g.vcf.gz",sample=SAMPLES)
    output:
        "data/work/IBD/sample_map.txt"
    run:
        with open(output[0],"w") as f:
            for sample in SAMPLES:
                f.write(f"{sample}\tdata/work/{sample}/gatk/common_snps.g.vcf.gz\n")

rule genomics_db_import:
    input:
        map="data/work/IBD/sample_map.txt",
        snps=f"data/work/common_snps/{id}/gnomad.exomes.v4.sites.{id}.common_biallelic_snps.simplified.vcf.gz"
    output:
        directory("data/work/IBD/snp_db")
    params:
        ref=config['referemce']['fasta']
    shell:
        """
        gatk GenomicsDBImport \
        --genomicsdb-workspace-path {output} \
        --sample-name-map {input.map} \
        -L {input.snps} \
        --merge-input-intervals false \
        --interval-padding 0
        """

rule genotype_gvcfs:
    input:
        db="data/work/IBD/snp_db",
        snps=f"data/work/common_snps/{id}/gnomad.exomes.v4.sites.{id}.common_biallelic_snps.simplified.vcf.gz"
    output:
        "data/work/IBD/joint.snp_genotypes.vcf.gz"
    params:
        ref=config['referemce']['fasta']
    shell:
        """
        gatk GenotypeGVCFs \
        -R {params.ref} \
        -V gendb://{input.db} \
        -L {input.snps} \
        --force-output-intervals {input.snps} \
        --include-non-variant-sites true \
        --call-genotypes true \
        -O {output}
        """

rule bcftools_vcf_norm:
    input:
        "data/work/IBD/joint.snp_genotypes.vcf.gz"
    output:
        "data/work/IBD/joint.snp_genotypes.norm.vcf.gz"
    shell:
        """
        bcftools norm -m-both -O z -o {output} {input}
        """

rule bcftools_vcf_clean:
    input:
        "data/work/IBD/joint.snp_genotypes.norm.vcf.gz"
    output:
        "data/work/IBD/joint.snp_genotypes.clean.vcf.gz"
    shell:
        """
        bcftools view -e 'ALT="*"' {input} |
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output}
        """

rule run_plink_genome:
    input:
        "data/work/IBD/joint.snp_genotypes.clean.vcf.gz"
    output:
        "data/work/IBD/plink.genome"
    params:
        out="data/work/IBD/plink"
    shell:
        """
        plink --vcf {input} --genome --out {params.out}
        """

rule report_ibd:
    input:
        "data/work/IBD/plink.genome"
    output:
        "data/work/IBD/ibd-related.txt"
    shell:
        """
        awk '$10 >= 0.1875 {{print $2, $4, $10}}' {input} > {output}
        """
##Make a better ibd report, with an R script, maybe even knit to html or something.

