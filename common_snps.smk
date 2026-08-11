# Canonical SNP resource for ibd/ibd2 (and metrics):
#   data/work/common_snps/{targets_key}/gnomad.exomes.common_biallelic_snps.{ref}.vcf.gz
#
# Set exactly one of these in config:
#
# common_snps — Pre-made common-SNP VCF (panel-wide or similar) that still needs
#   restriction to targets_intervals. Runs one GATK SelectVariants (-L intervals).
#
# gnomad_vcf — Broader gnomAD sites VCF; build common SNPs on the panel via
#   select -> simplify -> biallelic/AF filter -> final VCF at the path above.

rule common_biallelic_snps:
    input:
        expand("data/work/common_snps/{name}/gnomad.exomes.common_biallelic_snps.{ref}.vcf.gz",name=config["resources"]["targets_key"],ref=config["reference"]["key"])

if config["resources"]["common_snps"]:
    rule select_from_prefiltered:
        input:
            vcf=config["resources"]["common_snps"],
            intervals=config["resources"]["targets_intervals"]
        output:
            "data/work/common_snps/{name}/gnomad.exomes.common_biallelic_snps.{ref}.vcf.gz"
        shell:
            """
            gatk --java-options "-Xmx10g" SelectVariants \
            -V {input.vcf} \
            -L {input.intervals} \
            -O {output} \
            --lenient
            """

elif config["resources"]["gnomad_vcf"]:
    rule select_variants:
        input:
            vcf=config["resources"]["gnomad_vcf"],
            intervals=config["resources"]["targets_intervals"]
        output:
            temp("data/work/common_snps/{name}/selected.vcf")
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
            "data/work/common_snps/{name}/selected.vcf"
        output:
            temp("data/work/common_snps/{name}/simplified.vcf")
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
            "data/work/common_snps/{name}/simplified.vcf"
        output:
            temp("data/work/common_snps/{name}/simplified.vcf.idx")
        shell:
            """
            gatk --java-options "-Xmx10g" IndexFeatureFile -F {input}
            """

    rule filter_common_snps:
        input:
            vcf="data/work/common_snps/{name}/simplified.vcf",
            idx="data/work/common_snps/{name}/simplified.vcf.idx"
        output:
            "data/work/common_snps/{name}/gnomad.exomes.common_biallelic_snps.{ref}.vcf.gz"
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
else:
    raise ValueError("Configuration must specify either common_snps or gnomad_vcf")
