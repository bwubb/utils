rule common_biallelic_snps:
    input:
        expand("data/work/common_snps/{id}/gnomad.exomes.v4.sites.{id}.common_biallelic_snps.simplified.vcf.gz",id=config['resources']['targets_key'])

rule select_variants:
    input:
        vcf=config["resources"]["gnomad_vcf"],
        intervals=config["resources"]["intervals"]
    output:
        "data/work/common_snps/{id}/selected.vcf"
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants \
        -V {input.vcf} \
        -L {input.intervals} \
        -O {output} \
        --lenient
        """

rule simplify_vcf:
    input:
        "data/work/common_snps/{id}/selected.vcf"
    output:
        "data/work/common_snps/{id}/simplified.vcf"
    shell:
        """
        grep -v '#' {input} | grep -P "\\tPASS\\t" | \
        sed -e 's/PASS\\t.*AF=/PASS\\tAF=/g' -e 's/[;].*$//g' | \
        while read contig pos id ref alt qual filter info; do \
            printf "$contig\\t$pos\\t.\\t$ref\\t$alt\\t.\\t$filter\\t$info\\n"; \
        done > {output}
        """

rule index_simplified:
    input:
        "data/work/common_snps/{id}/simplified.vcf"
    output:
        "data/work/common_snps/{id}/simplified.vcf.idx"
    shell:
        """
        gatk --java-options "-Xmx10g" IndexFeatureFile -F {input}
        """

rule filter_common_snps:
    input:
        vcf="data/work/common_snps/{id}/simplified.vcf",
        idx="data/work/common_snps/{id}/simplified.vcf.idx"
    output:
        "data/work/common_snps/{id}/gnomad.exomes.v4.sites.{id}.common_biallelic_snps.simplified.vcf.gz"
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants \
        -V {input.vcf} \
        -select-type SNP \
        -restrict-alleles-to BIALLELIC \
        -select "AF > 0.01" \
        -O {output} \
        --lenient
        """
