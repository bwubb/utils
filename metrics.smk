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

header='SAMPLE BAIT_SET TOTAL_READS PCT_PF_READS PF_HQ_ALIGNED_READS PCT_READS_ALIGNED_IN_PAIRS'.split(' ')
if config['resources'].get('disambiguate',False):
    header+=['PCT_HUMAN','PCT_MOUSE','PCT_AMBIGUOUS']
header+=['PCT_DUP']
header+='PCT_SELECTED_BASES MEAN_TARGET_COVERAGE MEDIAN_TARGET_COVERAGE MAX_TARGET_COVERAGE PCT_USABLE_BASES_ON_TARGET ZERO_CVG_TARGETS_PCT PCT_TARGET_BASES_lt_20X PCT_TARGET_BASES_40X PCT_TARGET_BASES_100X'.split(' ')

localrules:summary,target_coverage_summary

rule standard_summary:
    input:
        expand('{project}.{targets}.{date}.metrics_summary.csv',project=config['project']['name'],targets=config['resources']['targets_key'],date=datetime.today().strftime('%Y%m%d')),
        #expand('{project}.{targets}.{date}.mean_target_coverage.csv',project=config['project']['name'],targets=config['resources']['targets_key'],date=datetime.today().strftime('%Y%m%d'))

###############################################

rule CollectAlignmentSummaryMetrics:
    input:
        sample_bam
    output:
        "metrics/{reference}/{sample}/alignment_summary.metrics"
    params:
        reference=config['reference']['fasta'],
        memory="10240m"
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar CollectAlignmentSummaryMetrics R={params.reference} I={input} O={output} VALIDATION_STRINGENCY=SILENT"

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
        #Default value = 0.05. but documentation states "If processing a small file, set the minimum percentage option (M) to 0.5, otherwise an error may occur."
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar CollectInsertSizeMetrics R={params.reference} I={input} O={output[0]} H={output[1]} M={params.MINIMUM_PCT} VALIDATION_STRINGENCY=SILENT"

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
    wildcard_constraints:
        target=config['resources']['targets_key']
    shell:
        "java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar CollectHsMetrics R={params.reference} I={input} O={output[0]} COVERAGE_CAP=1000 BAIT_INTERVALS={params.baits} TARGET_INTERVALS={params.targets} PER_TARGET_COVERAGE={output[1]} VALIDATION_STRINGENCY=SILENT"
#Message about new code syntax
#CollectHsMetrics -R /home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta -I bam_input/work/PP103-DZ2A/GRCh37/recal.bam -O metrics/S31285117/PP103-DZ2A/target.metrics -COVERAGE_CAP 1000 -BAIT_INTERVALS /home/bwubb/resources/Interval_files/SureSelect-Exon_v7.S31285117.Covered.picard.intervals -TARGET_INTERVALS /home/bwubb/resources/Interval_files/SureSelect-Exon_v7.S31285117.Covered.picard.intervals -PER_TARGET_COVERAGE metrics/S31285117/PP103-DZ2A/target_coverage.metrics -VALIDATION_STRINGENCY SILENT

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
        expand("metrics/{reference}/{sample}/alignment_summary.metrics",sample=SAMPLES,reference=config['reference']['key']),
        expand("metrics/{targets}/{sample}/target.metrics",sample=SAMPLES,targets=config['resources']['targets_key']),
        expand("metrics/{targets}/{sample}/insert_size.metrics",sample=SAMPLES,targets=config['resources']['targets_key']),
        expand("metrics/{reference}/{sample}/flagstat.metrics",sample=SAMPLES,reference=config['reference']['key'])

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

#This should be reduced to subset gene targets?
rule target_coverage_summary:
    input:
        expand("metrics/{targets}/{sample}/target_coverage.metrics",sample=SAMPLES,targets=config['resources']['targets_key'])
    output:
        expand('{project}.{targets}.{date}.mean_target_coverage.csv',project=config['project']['name'],targets=config['resources']['targets_key'],date=datetime.today().strftime('%Y%m%d'))
    run:
        target_coverage=defaultdict(lambda: defaultdict(float))
        sample_file=defaultdict(str)
        key_order=[]
        for sample in SAMPLES:
            sample_file[sample]=[x for x in input if sample in x][0]
            with open(sample_file[sample],'r') as cov_file:
                reader=csv.DictReader(cov_file,delimiter='\t')
                for row in reader:
                    key=tuple([row[x] for x in ['chrom','start','end','length','name']])
                    if key not in key_order:
                        key_order.append(key)
                    target_coverage[key][sample]=float(row['mean_coverage'])
        with open(output[0],'w') as outfile:
            writer=csv.writer(outfile,delimiter=',')
            writer.writerow(['Chr','Start','End','Length','Name']+SAMPLES)
            for key in key_order:
                row=list(key)+[target_coverage[key][sample] for sample in SAMPLES]
                writer.writerow(row)
