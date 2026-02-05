
rule mnp_regions:
    input:
        "data/final/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv"
    output:
        regions="data/work/mnp/chr{CHR}.regions.txt",
    run:
        import csv
        positions=set()
        with open(input.csv, newline="") as f:
            r=csv.DictReader(f)
            for row in r:
                chrom=row.get("Chr","").strip()
                pos_str=row.get("Start","").strip()
                if not chrom or not pos_str:
                    continue
                try:
                    pos=int(pos_str)
                    positions.add((chrom,pos))
                except ValueError:
                    pass
        with open(output.regions,'w') as f:
            for (chrom,pos) in sorted(positions):
                f.write(f"{chrom}\t{pos}\n")


rule bcftools_mnp_filter:
    input:
        vcf=config["mnp_vcf"],
        regions="data/work/mnp/chr{CHR}.regions.txt"
    output:
        bcf="data/work/mnp/chr{CHR}.mnp.bcf",
        csi="data/work/mnp/chr{CHR}.mnp.bcf.csi"
    params:
        id_file="data/work/mnp/chr{CHR}.id.txt"
    shell:
        """
        bcftools view -R {input.regions} -i 'ID=@{params.id_file}' {input.vcf} -W=csi -Ob -o {output}
        """

rule bcftools_mnp_sample_info:
    input:
        bcf="data/work/mnp/chr{CHR}.mnp.bcf"
    output:
        "data/work/mnp/chr{CHR}.mnp_sample_info.txt"
    shell:
        """
        bcftools query -f '[%ID %SAMPLE %GT %AD %DP %VAF\n]' {input} |
        egrep -v '.\/. | '0/0' > {output}
        """

rule test_mnp_sample_info:
    input:
        ""
        "data/work/mnp/chr{CHR}.mnp_sample_info.txt"
    output:
        "data/work/mnp/chr{CHR}.mnp_sample_info.PASS.txt",
        "data/work/mnp/chr{CHR}.mnp_sample_info.FAIL.txt"
    shell:
        "python test_mnp_sample_info.py {input} {output}"
    #This needs to be a kick ass script that reads the info and figures out 
    #If the samples have the variant pairs or not.
    #Im not sure how we should do this.