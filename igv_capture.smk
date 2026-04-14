import csv
import os

# Single-sample IGV screenshots. Input CSV columns (comma or tab):
#   Sample.ID, Chr, Start, End, Gene
# Optional config: delimiter (default ',')
#
# Outputs under: analysis/work/{sample}/somatic_variants/igv_png/
#   {mut}.bwa.bam, {mut}.mutect2.bam, .bat, .png
# bam.table: sample<TAB>path per line
# Set IGV jar path in rule run_bat if needed.

MUT_LOOKUP={}

def _input_delim():
    return config.get('delimiter', ',')

def get_position(wildcards):
    entry=MUT_LOOKUP[f'{wildcards.sample}_{wildcards.mut}']
    return entry['range']

def igv_snapshot_bam(wildcards):
    base=f"analysis/work/{wildcards.sample}/somatic_variants/igv_png/{wildcards.mut}"
    return f"{base}.mutect2.bam"

def igv_snapshot_bai(wildcards):
    base=f"analysis/work/{wildcards.sample}/somatic_variants/igv_png/{wildcards.mut}"
    return f"{base}.mutect2.bam.bai"

MUTS={}
RANGES={}
with open(config['input'],'r',newline='') as file:
    reader=csv.DictReader(file,delimiter=_input_delim())
    for n,row in enumerate(reader):
        k=f'n{n}'
        sample=row['Sample.ID'].strip()
        chrom=row['Chr'].strip()
        start=int(row['Start'])
        end=int(row['End'])
        gene=row['Gene'].strip()
        os.makedirs(f'analysis/work/{sample}/somatic_variants/igv_png',exist_ok=True)
        # unique per row (sample + locus + gene)
        mut_name=f'{sample}.{chrom}-{start}-{end}.{gene}'
        # samtools / GATK region (contig:from-to)
        region=f'{chrom}:{start}-{end}'
        # IGV goto (same interval; extend slightly if you prefer)
        igv_range=region
        MUTS[k]={'name':mut_name,'range':igv_range,'sample':sample}
        bwa_path=f'analysis/work/{sample}/somatic_variants/igv_png/{mut_name}.bwa.bam'
        RANGES[bwa_path]=region
        MUT_LOOKUP[f'{sample}_{mut_name}']=MUTS[k]

TARGET_FILES=[]
SCREENSHOT_FILES=[]
for k,v in MUTS.items():
    mut_name=v['name']
    sample=v['sample']
    TARGET_FILES.append(f'analysis/work/{sample}/somatic_variants/igv_png/{mut_name}.bwa.bam')
    TARGET_FILES.append(f'analysis/work/{sample}/somatic_variants/igv_png/{mut_name}.mutect2.bam')
    SCREENSHOT_FILES.append(f'analysis/work/{sample}/somatic_variants/igv_png/{mut_name}.bat')
    SCREENSHOT_FILES.append(f'analysis/work/{sample}/somatic_variants/igv_png/{sample}_{mut_name}.png')

BAMS={}
with open('bam.table','r') as file:
    for line in file:
        sample,bam=line.rstrip().split('\t')
        BAMS[sample]=bam

wildcard_constraints:
    mut='[^/]+'

rule all:
    input:
        TARGET_FILES+SCREENSHOT_FILES

rule samtools_bamout:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        bam="analysis/work/{sample}/somatic_variants/igv_png/{mut}.bwa.bam",
        bai="analysis/work/{sample}/somatic_variants/igv_png/{mut}.bwa.bam.bai",
    params:
        range=lambda wildcards: RANGES[f"analysis/work/{wildcards.sample}/somatic_variants/igv_png/{wildcards.mut}.bwa.bam"]
    shell:
        """
        samtools view -O BAM -o {output.bam} {input.bam} {params.range}
        samtools index {output.bam}
        """

rule mutect2_bamout:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        bam="analysis/work/{sample}/somatic_variants/igv_png/{mut}.mutect2.bam",
        bai="analysis/work/{sample}/somatic_variants/igv_png/{mut}.mutect2.bam.bai",
        vcf="analysis/work/{sample}/somatic_variants/igv_png/{mut}.mutect2.vcf"
    params:
        ref=config['reference']['fasta'],
        sample_id=lambda wildcards: wildcards.sample,
        range=lambda wildcards: RANGES[f"analysis/work/{wildcards.sample}/somatic_variants/igv_png/{wildcards.mut}.bwa.bam"]
    log:
        "analysis/work/{sample}/somatic_variants/igv_png/{mut}.mutect2.log"
    shell:
        """
        gatk Mutect2 -R {params.ref} \
        -I {input.bam} \
        -tumor {params.sample_id} \
        -L {params.range} \
        -ip 1000 -bamout {output.bam} -O {output.vcf}

        samtools index {output.bam}
        """

rule write_bat:
    input:
        bam=igv_snapshot_bam,
        bai=igv_snapshot_bai,
    output:
        bat="analysis/work/{sample}/somatic_variants/igv_png/{mut}.bat"
    params:
        snapshot_dir=lambda wildcards: f"analysis/work/{wildcards.sample}/somatic_variants/igv_png",
        snapshot_name=lambda wildcards: f"{wildcards.sample}_{wildcards.mut}",
        position=get_position,
    run:
        with open(output['bat'],'w') as file:
            file.write(f"new\n")
            file.write(f"genome {config['reference']['key']}\n")
            file.write(f"snapshotDirectory {params.snapshot_dir}/\n")
            file.write(f"load {input['bam']}\n")
            file.write(f"maxPanelHeight 2000\n")
            file.write(f"goto {params.position}\n")
            file.write(f"group strand\n")
            file.write(f"snapshot {params.snapshot_name}\n")
            file.write(f"exit\n")

rule run_bat:
    input:
        bat="analysis/work/{sample}/somatic_variants/igv_png/{mut}.bat",
        bam=igv_snapshot_bam,
        bai=igv_snapshot_bai,
    output:
        png="analysis/work/{sample}/somatic_variants/igv_png/{sample}_{mut}.png"
    shell:
        "xvfb-run --auto-servernum --server-num=1 --server-args='-screen 0, 1900x1200x24' java -Xmx4g -jar /usr/local/software/IGV_2.4.4/igv.jar -b {input.bat}"
