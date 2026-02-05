import os
import sys
import yaml

# Create log directories
os.makedirs('logs/cluster/aln',exist_ok=True)
os.makedirs('logs/cluster/metrics',exist_ok=True)


### FUNCTIONS ###

def map_input(wildcards):
    inputs=[]
    for RUN,_run in list(FILES[wildcards.sample].items()):
        run,lane,index=_run['PU'].split('-',2)#this causes annoyances when - is used in run id. Not fixed by rsplit.
        inputs.append(f'bam_input/work/{wildcards.sample}/{wildcards.reference}/{run}/{lane}/{index}/5.markdup.bam')
    assert len(inputs)>0
    return sorted(inputs)

def get_fastqs(wildcards):
    key=f'{wildcards.run}-{wildcards.lane}-{wildcards.index}'
    entry=FILES[wildcards.sample][key]
    r1,r2=entry['files'][0],entry['files'][1]
    def p(x):
        return x if x.startswith('FASTQ') or x.startswith('/') else 'FASTQ/'+x
    return {'R1':p(r1),'R2':p(r2)}

### ### PYTHON ### ###

with open(config['project']['fastq_config']) as file:
    FILES=yaml.load(file,Loader=yaml.BaseLoader)
    for sample in list(FILES.keys()):
        val=FILES[sample]
        flist=sorted(val['files'])
        r1_by_base={}
        r2_by_base={}
        for f in flist:
            base=os.path.basename(f)
            if '_R1' in base:
                base_stem=base.replace('_R1.fastq.gz','').replace('_R1.fq.gz','')
                r1_by_base[base_stem]=f
            elif '_R2' in base:
                base_stem=base.replace('_R2.fastq.gz','').replace('_R2.fq.gz','')
                r2_by_base[base_stem]=f
        paired=[]
        for base_stem in sorted(set(r1_by_base)|set(r2_by_base)):
            r1=r1_by_base.get(base_stem)
            r2=r2_by_base.get(base_stem)
            if r1 is None:
                print(f"WARNING: unpaired R2 (no R1): {r2}",file=sys.stderr)
                continue
            if r2 is None:
                print(f"WARNING: unpaired R1 (no R2): {r1}",file=sys.stderr)
                continue
            paired.append((base_stem,r1,r2))
        if not paired:
            print(f"WARNING: no paired FASTQs for sample {sample}, skipping",file=sys.stderr)
            continue
        new_val={}
        for base_stem,r1,r2 in paired:
            parts=base_stem.split('_')
            run,lane,index=(parts[-3],parts[-2],parts[-1]) if len(parts)>=3 else ('run','0','0')
            key=f'{run}-{lane}-{index}'
            new_val[key]={'PU':f'{run}-{lane}-{index}','files':[r1,r2]}
        FILES[sample]=new_val

with open(config['project']['sample_list']) as file:
    SAMPLES=file.read().splitlines()
for sample in SAMPLES:
    if sample not in FILES:
        print(f"WARNING: sample {sample} in sample.list not in fastq config, skipping",file=sys.stderr)
SAMPLES=[s for s in SAMPLES if s in FILES]
for sample in SAMPLES:
    os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

### ### ### RULES ### ### ###

localrules:aln_all,bam_table,write_bam_table

rule aln_all:
    input:
        expand("bam_input/final/{sample}/{sample}.{reference}.bam",sample=SAMPLES,reference=config['reference']['key'])

rule bam_table:
    input:
        "bam.table"

rule bwa_mem:
    input:
        unpack(get_fastqs)
    output:
        temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"),
    params:
        fasta=config['reference']['fasta']
    threads:
        4
    shell:
        #fastq config has no path.
        """
        bwa mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools view -bS -o {output}
        """

rule samtools_readgroup:
    input:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"
    output:
        rg=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/2.readgroup.bam"),
        fm=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/3.fixmate.bam"),
        qs=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/4.qsort.bam"),
    params:
        LB=config['resources']['library_key']
    threads:
        4
    shell:
        """
        samtools addreplacerg -@ {threads} -r 'ID:{wildcards.run}.{wildcards.lane}' -r 'PU:{wildcards.run}.{wildcards.lane}.{wildcards.index}' -r 'PL:illumina' -r 'LB:{params.LB}' -r 'SM:{wildcards.sample}' -o {output.rg} {input}
        samtools fixmate -m -@ {threads} {output.rg} {output.fm}
        samtools sort -@ {threads} -o {output.qs} {output.fm}
        """

rule samtools_markdup:
    input:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/4.qsort.bam"
    output:
        temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/5.markdup.bam")
    params:
        stats="bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/5.stats.txt"
    threads:
        4
    shell:
        """
        samtools markdup -s -f {params.stats} -@ {threads} {input} {output}
        """

rule input_ready:
    input:
        map_input
    output:
        temp("bam_input/work/{sample}/{reference}/input.bam")
    threads:
        4
    shell:
        """
        samtools merge -f -@ {threads} {output} {input}
        samtools index {output}
        """

rule ValidateSamFile:
    input:
        "bam_input/work/{sample}/{reference}/input.bam"
    output:
        "metrics/{reference}/{sample}/validation_data.table"
    params:
        memory="10240m"
    shell:
        """
        set +e
        java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar ValidateSamFile I={input} O={output} MODE=SUMMARY
        exitcode=$?
        if [ $exitcode -ne 0 ]; then exit 1; fi
        exit 0
        """
        #Is this also in GATK? Then I could drop picard from this pipeline.

#This is very strict
#MATE_CIGAR_STRING_INVALID_PRESENCE
#May need case when to allow some errors.
#rule validation_pass:
#    input:
#        "bam_input/work/{sample}/{reference}/validation_data.table"
#    output:
#        "metrics/{reference}/{sample}/validation_data.table"
#    shell:
#        """
#        set +H
#        if egrep -q 'No errors found' {input[0]}; then
#            cp {input[0]} {output[0]}
#        else
#            egrep '^ERROR' {input[0]}
#            exit 1
#        fi
#        """

rule ready_bam:
    input:
        bam="bam_input/work/{sample}/{reference}/input.bam",
        table="metrics/{reference}/{sample}/validation_data.table"
    output:
        "bam_input/final/{sample}/{sample}.{reference}.bam"
    shell:
        """
        rsync -v {input.bam} {output}
        samtools index {output}
        """

rule write_bam_table:
    input:
        expand("bam_input/final/{sample}/{sample}.{reference}.bam",sample=SAMPLES,reference=config['reference']['key'])
    output:
        "bam.table"
    params:
        ref=config['reference']['key']
    run:
        with open(output[0],'w') as file:
            for sample in SAMPLES:
                file.write(f"{sample}\tbam_input/final/{sample}/{sample}.{params.ref}.bam\n")
