import csv
import os
from collections import defaultdict

def mkdir_p(path):
    os.makedirs(path,exist_ok=True)

def get_position(wildcards):
    position=MUTS[f'{wildcards.tumor}_{wildcards.mut}']
    chr,start=position.split(':')
    return f'{chr}:{int(start)-100}-{int(start)+100}'

def paired_bams(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

MUTS={}
RANGES={}
with open(config['input'],'r') as file:
    #dict['GENE_pA123B']=1:23345
    reader=csv.DictReader(file,delimiter='\t')
    for n,row in enumerate(reader):
        k=f'n{n}'
        if 'Tumor.ID' in row.keys():
            samples=[row['Tumor.ID'],row['Normal.ID']]
            os.makedirs(f'analysis/work/{row["Tumor.ID"]}/somatic_variants/igv_png',exist_ok=True)
            os.makedirs(f'analysis/work/{row["Normal.ID"]}/somatic_variants/igv_png',exist_ok=True)
        else:
            samples=[row['Sample.ID']]
            os.makedirs(f'analysis/work/{row["Sample.ID"]}/somatic_variants/igv_png',exist_ok=True)
        chr=row['Chr']
        start=row['Start']
        gene=row['Gene']
        if row.get('HGVSp',f'variant{n}')=='.':
            hgvs='splice'
        else:
            hgvs=row.get('HGVSp',f'variant{n}').replace('p.','').replace('?','')
        MUTS[k]={'name':f'chr{chr}-{start}.{gene}.{hgvs}','range':f'{chr}:{int(start)-100}-{int(start)+100}','samples':samples}
        for sample in samples:
            #MUTS[f'analysis/work/{sample}/chr{chr}-{start}.{gene}.{hgvs}.bwa.bam']
            RANGES[f'analysis/work/{sample}/chr{chr}-{start}.{gene}.{hgvs}.bwa.bam']=f'{chr}:{int(start)-100}-{int(start)+100}'

TARGET_FILES=[]
for k,v in MUTS.items():
    for sample in v['samples']:
        TARGET_FILES.append(f'analysis/work/{sample}/{v["name"]}.bwa.bam')
        TARGET_FILES.append(f'analysis/work/{sample}/{v["name"]}.mutect2.bam')

BAMS={}
with open('bam.table','r') as file:
    for line in file:
        sample,bam=line.rstrip().split('\t')
        BAMS[sample]=bam

PAIRS={}
with open('pair.table','r') as file:
    for line in file:
        tumor,normal=line.rstrip().split('\t')
        PAIRS[tumor]=normal

wildcard_constraints:
    name='chr\d+-\d+.[a-zA-Z0-9_.]+.[a-zA-Z0-9_.]+'

rule all:
    input:
        TARGET_FILES

rule samtools_bamout:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "analysis/work/{sample}/{mut}.bwa.bam"
    params:
        range=lambda wildcards: RANGES[f"analysis/work/{wildcards.sample}/{wildcards.mut}.bwa.bam"]
    shell:
        """
        samtools view -O BAM -o {output} {input.bam} {params.range}
        samtools index {output}
        """

rule mutect2_bamout:
    input:
        unpack(paired_bams)
    output:
        bam="analysis/work/{tumor}/{mut}.mutect2.bam",
        vcf="analysis/work/{tumor}/{mut}.mutect2.vcf"
    params:
        ref=config['reference']['fasta'],
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        range=lambda wildcards: RANGES[f"analysis/work/{wildcards.tumor}/{wildcards.mut}.bwa.bam"]
    log:
        "analysis/work/{tumor}/{mut}.mutect2.log"
    shell:
        """
        gatk Mutect2 -R {params.ref} \
        -I {input.tumor} -I {input.normal} \
        -tumor {params.tumor} -normal {params.normal} \
        -L {params.range} \
        -ip 1000 -bamout {output.bam} -O {output.vcf}

        samtools index {output.bam}
        """

#gatk Mutect2 -R {params.ref} -I {input.tumor} -I {input.normal} \
#-tumor {params.tumor} -normal {params.normal} -L {params.intervals} \
#-ip 1000 -bamout {output.bam} -O {output.vcf}
#bam="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.bam",
#vcf="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2_activeregion.vcf.gz"


#rule annotate_activeregion_vcf:
#    input:
#        "analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2_activeregion.vcf.gz"
#    output:
#        "analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.vcf.gz"
#    shell:
#        """
#        bcftools view
#        """

rule write_bat:
    input:
        bam="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.bam"
    output:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.bat"
    params:
        dir="analysis/work/{tumor}/somatic_variants/igv_png/",
        name="{tumor}_{mut}",
        position=get_position
    run:
        with open(output['bat'],'w') as file:
            file.write(f"new\n")
            file.write(f"genome {config['reference']['key']}\n")
            file.write(f"snapshotDirectory {params['dir']}\n")
            file.write(f"load {input['bam']}\n")
            file.write(f"maxPanelHeight 2000\n")
            file.write(f"goto {params['position']}\n")
            file.write(f"group strand\n")
            file.write(f"snapshot {params['name']}\n")
            file.write(f"exit\n")

rule run_bat:
    input:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.bat"
    output:
        png="analysis/work/{tumor}/somatic_variants/igv_png/{tumor}_{mut}.png"
    shell:
        "xvfb-run --auto-servernum --server-num=1 --server-args='-screen 0, 1900x1200x24' java -Xmx4g -jar /usr/local/software/IGV_2.4.4/igv.jar -b {input}"
