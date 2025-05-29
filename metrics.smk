import os
import csv
from datetime import datetime
from collections import defaultdict

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()

with open(config.get('project',{}).get('bam_table','bams.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def sample_bam(wildcards):
    return BAMS[wildcards.sample]

# Define headers in a cleaner way
base_header=['SAMPLE','BAIT_SET','TOTAL_READS','PCT_PF_READS','PF_HQ_ALIGNED_READS','PCT_READS_ALIGNED_IN_PAIRS']
if config['resources'].get('disambiguate',False):
    base_header.extend(['PCT_HUMAN','PCT_MOUSE','PCT_AMBIGUOUS'])
base_header.extend(['PCT_DUP'])
base_header.extend(['PCT_SELECTED_BASES','MEAN_TARGET_COVERAGE','MEDIAN_TARGET_COVERAGE','MAX_TARGET_COVERAGE',
                   'PCT_USABLE_BASES_ON_TARGET','ZERO_CVG_TARGETS_PCT','PCT_TARGET_BASES_lt_20X',
                   'PCT_TARGET_BASES_40X','PCT_TARGET_BASES_100X'])

localrules: metrics_summary,target_coverage_summary

rule standard_summary:
    input:
        expand('{project}.{targets}.{date}.metrics_summary.csv',
               project=config['project']['name'],
               targets=config['resources']['targets_key'],
               date=datetime.today().strftime('%Y%m%d'))

rule CollectAlignmentSummaryMetrics:
    input:
        sample_bam
    output:
        "metrics/{reference}/{sample}/alignment_summary.metrics"
    params:
        reference=config['reference']['fasta'],
        memory="10240m"
    shell:
        "gatk --java-options -Xmx{params.memory} CollectAlignmentSummaryMetrics -R {params.reference} -I {input} -O {output} --VALIDATION_STRINGENCY SILENT"

rule CollectInsertSizeMetrics:
    input:
        sample_bam
    output:
        "metrics/{targets}/{sample}/insert_size.metrics",
        "metrics/{targets}/{sample}/insert_size_histogram.pdf"
    params:
        reference=config['reference']['fasta'],
        memory="10240m",
        MINIMUM_PCT=0.05
    shell:
        "gatk --java-options -Xmx{params.memory} CollectInsertSizeMetrics -R {params.reference} -I {input} -O {output[0]} -H {output[1]} -M {params.MINIMUM_PCT} --VALIDATION_STRINGENCY SILENT"

rule CollectHsMetrics:
    input:
        sample_bam
    output:
        "metrics/{targets}/{sample}/target.metrics",
        "metrics/{targets}/{sample}/target_coverage.metrics"
    params:
        reference=config['reference']['fasta'],
        baits=config['resources']['picard_intervals'],
        targets=config['resources']['picard_intervals'],
        memory="10240m"
    shell:
        "gatk --java-options -Xmx{params.memory} CollectHsMetrics -R {params.reference} -I {input} -O {output[0]} --COVERAGE_CAP 1000 --BAIT_INTERVALS {params.baits} --TARGET_INTERVALS {params.targets} --PER_TARGET_COVERAGE {output[1]} --VALIDATION_STRINGENCY SILENT"

rule samtools_flagstat:
    input:
        sample_bam
    output:
        "metrics/{reference}/{sample}/flagstat.metrics"
    shell:
        "samtools flagstat {input} > {output}"

rule samtools_readlength:
    input:
        sample_bam
    output:
        "metrics/{reference}/{sample}/readlength.metrics"
    shell:
        "samtools view -F 4 {input} | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . \"\n\"' | sort -n | uniq -c > {output}"

rule metrics_summary:
    input:
        alignment=expand("metrics/{reference}/{sample}/alignment_summary.metrics",
                        sample=SAMPLES,
                        reference=config['reference']['key']),
        target=expand("metrics/{targets}/{sample}/target.metrics",
                     sample=SAMPLES,
                     targets=config['resources']['targets_key']),
        insert_size=expand("metrics/{targets}/{sample}/insert_size.metrics",
                     sample=SAMPLES,
                     targets=config['resources']['targets_key']),
        flagstat=expand("metrics/{reference}/{sample}/flagstat.metrics",
                       sample=SAMPLES,
                       reference=config['reference']['key'])
    output:
        csv='{project}.{date}.metrics_summary.csv'
    params:
        PDX=str('disambiguate' in config),
        sample_list=config['project']['sample_list'],
        ref=config['reference']['key'],
        targets=config['resources']['targets_key']
    shell:
        """
        python metrics_summary.py -I {params.sample_list} -O {output} -L {params.targets} -R {params.ref} --PDX {params.PDX}
        """

rule target_coverage_summary:
    input:
        expand("metrics/{targets}/{sample}/target_coverage.metrics",
               sample=SAMPLES,
               targets=config['resources']['targets_key'])
    output:
        expand('{project}.{targets}.{date}.mean_target_coverage.csv',
               project=config['project']['name'],
               targets=config['resources']['targets_key'],
               date=datetime.today().strftime('%Y%m%d'))
    params:
        sample_list=config['project']['sample_list'],
        targets=config['resources']['targets_key']
    shell:
        """
        python target_coverage_summary.py -I {params.sample_list} -O {output} -L {params.targets}
        """
