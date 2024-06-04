import csv
import os
from collections import defaultdict

#Read tumor/normals names/bams
#Read mutation file
#Make intervals
#Make bed +500 each side
#IGV thing to make bat, -nosnap, Can I make my own with better settings???
#Run xvfb with modified settings
#???
#Profit


def get_position(wildcards):
    position=MUTS[wildcards.mut]
    chr,start=position.split(':')
    return f'{chr}:{int(start)-100}-{int(start)+100}'

def get_window(wildcards):
    position=MUTS[wildcards.mut]
    chr,start=position.split(':')
    return f'{chr}:{int(start)-1000}-{int(start)+1000}'

def get_bam(wildcards):
    return BAMS[wildcards.sample]

#TumorID', 'MayoBrca2Br34
MUTS={}
SAMPLE_MUT=defaultdict(list)
with open('','r') as file:
    #dict['GENE_pA123B']=1:23345
    reader=csv.DictReader(file,delimiter=',')
    for row in reader:
        #if row['AAChange.refGene']=='.':
        #    X=row['NTChange.refGene'].replace('.','')
        #    #If not AAChange.refGene use NTChange.refGene
        #    #Need to check for special characters
        #else:
        #    X=row['AAChange.refGene'].replace('.','')
        #key=f"{row['Gene.refGene']}_{X}"
        #I have no standardized naming convention.
        key=row['SpliceAI']

        MUTS[key]=f"{row['Chr']}:{row['Start']}"
        SAMPLE_MUT[row['SampleID']].append(key)
        os.makedirs(f"analysis/work/{row['SampleID']}/germline_variants/igv_png",exist_ok=True)
        with open(f"analysis/work/{row['SampleID']}/germline_variants/igv_png/{key}.bed",'w') as bed_file:
            bed_row=[row['Chr']]
            bed_row.append(str(int(row['Start'])-100))
            bed_row.append(str(int(row['Start'])+100))
            bed_row.append(key)
            out="\t".join(bed_row)
            bed_file.write(f'{out}\n')

TARGET_FILES=[]
TARGET_BAMS=[]
for sample,muts in SAMPLE_MUT.items():
    for mut in muts:
        TARGET_FILES.append(f'analysis/work/{sample}/germline_variants/igv_png/{sample}_{mut}.png')
        TARGET_BAMS.append(f"analysis/work/{sample}/germline_variants/igv_png/{sample}_{mut}.gatk.bam")
        TARGET_BAMS.append(f"analysis/work/{sample}/germline_variants/igv_png/{sample}_{mut}.bwa.bam")



BAMS={}
with open('bam.table','r') as file:
    for line in file:
        sample,bam=line.rstrip().split('\t')
        BAMS[sample]=bam

rule variant_bams:
    input:
        TARGET_BAMS

rule all:
    input:
        TARGET_FILES

#look at --assembly-region-out and --activity-profile-out

rule GATK_bamout:
    input:
        get_bam
    output:
        bam="analysis/work/{sample}/germline_variants/igv_png/{sample}_{mut}.gatk.bam",
        vcf="analysis/work/{sample}/germline_variants/igv_png/{mut}.gatk_activeregion.vcf.gz"
    params:
        ref='/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta',
        intervals=lambda wildcards: MUTS[wildcards.mut]
    log:
        "analysis/work/{sample}/germline_variants/igv_png/{mut}.gatk.log"
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller \
        -R {params.ref} -I {input} \
        -L {params.intervals} \
        -ip 1000 -bamout {output.bam} -O {output.vcf}

        samtools index {output.bam}
        """
#--assembly-region-out
#	null 	Output the assembly region to this IGV formatted file
#--base-quality-score-threshold
#	18 	Base qualities below this threshold will be reduced to the minimum (6)
#--disable-tool-default-read-filters
#--force-active

#rule annotate_activeregion_vcf:
#    input:
#        "analysis/work/{sample}/germline_variants/igv_png/{mut}.mutect2_activeregion.vcf.gz"
#    output:
#        "analysis/work/{sample}/germline_variants/igv_png/{mut}.mutect2.vcf.gz"
#    shell:
#        """
#        bcftools view
#        """
rule BWA_bamout:
    input:
        get_bam
    output:
        bam="analysis/work/{sample}/germline_variants/igv_png/{sample}_{mut}.bwa.bam"
    params:
        ref='/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta',
        window=get_window
    shell:
        """
        samtools view -b -o {output} {input} {params.window}
        samtools index {output.bam}
        """

rule write_bat:
    input:
        bam1="analysis/work/{sample}/germline_variants/igv_png/{mut}.gatk.bam",
        bam2="analysis/work/{sample}/germline_variants/igv_png/{mut}.bwa.bam"
    output:
        bat="analysis/work/{sample}/germline_variants/igv_png/{mut}.bat"
    params:
        dir="analysis/work/{sample}/germline_variants/igv_png/",
        name="{sample}_{mut}",
        position=get_position
    run:
        with open(output['bat'],'w') as file:
            file.write(f"new\n")
            file.write(f"genome b37\n")
            file.write(f"snapshotDirectory {params['dir']}\n")
            file.write(f"load {input['bam1']}\n")
            file.write(f"load {input['bam2']}\n")
            file.write(f"maxPanelHeight 2000\n")
            file.write(f"goto {params['position']}\n")
            file.write(f"group strand\n")
            file.write(f"snapshot {params['name']}\n")
            file.write(f"exit\n")

rule run_bat:
    input:
        bat="analysis/work/{sample}/germline_variants/igv_png/{mut}.bat"
    output:
        bat="analysis/work/{sample}/germline_variants/igv_png/{sample}_{mut}.png"
    shell:
        "xvfb-run --auto-servernum --server-num=1 --server-args='-screen 0, 1900x1200x24' java -Xmx10240m -jar /usr/local/software/IGV_2.4.4/igv.jar -b {input}"
