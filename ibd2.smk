import os
import yaml

include: "common_snps.smk"
include: "metrics.smk"

os.makedirs("logs/cluster/ibd2",exist_ok=True)

with open(config.get('project',{}).get('sample_list','sample.list'),'r') as i:
    SAMPLES=i.read().splitlines()

with open(config.get('project',{}).get('bam_table','bams.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

name=config['resources']['targets_key']
ref=config['reference']['key']

# Define chromosomes
CHROMOSOMES=[f'chr{i}' for i in range(1,23)]+['chrX','chrY']

def get_gvcf_variants():
    """Return list of gVCF files for GATK CombineGVCFs --variant parameter"""
    return [f"data/work/{sample}/gatk/common_snps.g.vcf.gz" for sample in SAMPLES]

localrules: make_sample_map,combine_gvcfs_all,genotype_gvcfs_all

rule collect_ibd:
    input:
        "data/work/IBD/ibd-related.txt"

rule gatk_alleles:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample],
        intervals="data/work/common_snps/{name}/snps.intervals",
        snp="data/work/common_snps/{name}/gnomad.exomes.common_biallelic_snps.{ref}.vcf.gz"
    output:
        "data/work/gatk/{sample}/common_snps.g.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        gatk --java-options '-Xmx5g' HaplotypeCaller \
        -R {params.ref} \
        -I {input.bam} \
        -O {output} \
        -L {input.intervals} \
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
        "data/work/IBD2/sample_map.txt"
    run:
        with open(output[0],"w") as f:
            for sample in SAMPLES:
                f.write(f"{sample}\tdata/work/{sample}/gatk/common_snps.g.vcf.gz\n")

#I think the query makes the small intervals file, then the variant arguments are the gvcfs
#at no point have we actually subset anything to chromosome level.
rule combine_gvcfs_chr:
    input:
        gvcfs=expand("data/work/{sample}/gatk/common_snps.g.vcf.gz",sample=SAMPLES),
        snps=f"data/work/common_snps/{name}/gnomad.exomes.common_biallelic_snps.{ref}.vcf.gz"
    output:
        "data/work/IBD2/combined.{chr}.g.vcf.gz"
    params:
        ref=config['reference']['fasta'],
        chr="{chr}"
    shell:
        """
        # Extract chromosome intervals from SNP file
        bcftools query -f '%CHROM:%POS-%POS\n' {input.snps} | grep ^{params.chr}: > {params.chr}.intervals
        
        # Generate variant arguments
        variant_args=""
        for gvcf in {input.gvcfs}; do
            variant_args="$variant_args --variant $gvcf"
        done
        
        # Combine GVCFs for this chromosome
        gatk --java-options '-Xmx16g' CombineGVCFs \
        --reference {params.ref} \
        $variant_args \
        --intervals {params.chr}.intervals \
        --output {output}
        
        # Cleanup
        rm {params.chr}.intervals
        """

rule combine_gvcfs_all:
    input:
        expand("data/work/IBD2/combined.{chr}.g.vcf.gz",chr=CHROMOSOMES)

rule genotype_gvcfs_chr:
    input:
        "data/work/IBD2/combined.{chr}.g.vcf.gz"
    output:
        "data/work/IBD2/joint.{chr}.snp_genotypes.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        gatk --java-options '-Xmx16g' GenotypeGVCFs \
        -R {params.ref} \
        -V {input} \
        -O {output}
        """

rule genotype_gvcfs_all:
    input:
        expand("data/work/IBD2/joint.{chr}.snp_genotypes.vcf.gz",chr=CHROMOSOMES)

rule merge_chromosomes:
    input:
        expand("data/work/IBD2/joint.{chr}.snp_genotypes.vcf.gz",chr=CHROMOSOMES)
    output:
        "data/work/IBD2/joint.snp_genotypes.vcf.gz"
    shell:
        """
        # Merge all chromosomes
        bcftools concat -a -D -O z -o {output} {input}
        bcftools index {output}
        """

rule bcftools_vcf_norm:
    input:
        "data/work/IBD2/joint.snp_genotypes.vcf.gz"
    output:
        "data/work/IBD2/joint.snp_genotypes.norm.vcf.gz"
    shell:
        """
        bcftools norm -m-both -O z -o {output} {input}
        """

rule bcftools_vcf_clean:
    input:
        "data/work/IBD2/joint.snp_genotypes.norm.vcf.gz"
    output:
        "data/work/IBD2/joint.snp_genotypes.clean.vcf.gz"
    shell:
        """
        bcftools view -e 'ALT="*"' {input} |
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output}
        """

rule run_plink_genome:
    input:
        "data/work/IBD2/joint.snp_genotypes.clean.vcf.gz"
    output:
        "data/work/IBD2/plink.genome"
    params:
        out="data/work/IBD2/plink"
    shell:
        """
        plink --vcf {input} --genome --out {params.out}
        """

rule report_ibd:
    input:
        "data/work/IBD2/plink.genome"
    output:
        "data/work/IBD2/ibd-related.txt",
        "data/work/IBD2/ibd_report.html"
    shell:
        """
        awk '$10 >= 0.1875 {{print $2, $4, $10}}' {input} > {output[0]}
        Rscript -e 'library(rmarkdown); render("ibd_report.Rmd", output_file="{output[1]}", params=list(genome_file="{input}"))'
        """
