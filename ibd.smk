

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()

with open(config.get('project',{}).get('bam_table','bams.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())
#I need to standardize name the gnomad stuff.

rule collect_ibd:
    input:
        "data/work/IBD/ibd-related.txt"

rule gatk_alleles:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
        #/home/bwubb/resources/Vcf_files/gnomad.exomes.liftover_grch38.S31285117.common_biallelic_snps.simplified.vcf.gz
    output:
        "data/work/{lib}/{sample}/gatk/haplotype.alleles.vcf.gz"
    params:
        ref=config["reference"]["fasta"],
        alleles=config["resources"]["common_snps"]
    shell:
        """
        gatk --java-options '-Xmx10240m' HaplotypeCaller \
        -I {input.bam} \
        -O {output} \
        -R {params.ref} \
        --alleles {params.alleles} \
        """

rule make_input_list:
    input:
        expand("data/work/{lib}/{sample}/gatk/haplotype.alleles.vcf.gz",sample=SAMPLES,lib=config["resources"]["targets_key"])
    output:
        "data/work/IBD/input.list"
    run:
        with open(f"{output[0]}","w") as file:
            for i in input:
                file.write(f"{i}\n")

#trim down to alleles
#remove *
#merge


rule bcftools_vcf_merge:
    input:
        "data/work/IBD/input.list"
    output:
        "data/work/IBD/merged.haplotype.alleles.1.vcf.gz"
    params:
        alleles=config["resources"]["common_snps"]
    shell:
        """
        bcftools merge \
        -l {input} \
        -R {params.alleles} \
        -O z \
        -o {output} \
        """

rule run_plink_genome:
    input:
        "data/work/IBD/merged.haplotype.alleles.1.vcf.gz"
    output:
        "data/work/IBD/plink.genome"
    shell:
        """
        plink --vcf {input} --genome --out {output}
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
