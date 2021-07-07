import glob
import os
import csv
import argparse
from collections import defaultdict

def flagstat_metrics(FileO):
    TOTAL=int(FileO.readline().rstrip().split(' ')[0])
    for line in FileO:
        if "duplicates" in line:
            DUP=int(line.split(' ')[0])
            break
    try:
        return {'PCT_DUP':str(float(DUP/TOTAL))}
    except ZeroDivisionError as e:
        return {'PCT_DUP':'.'}

def alignment_metrics(FileO):
    for line in FileO:
        if line.startswith('CATEGORY'):
            break
    reader=csv.DictReader(FileO,delimiter='\t',fieldnames=line.rstrip().split('\t'))
    FIRST_OF_PAIR=reader.__next__()
    SECOND_OF_PAIR=reader.__next__()
    PAIR=reader.__next__()
    return PAIR

def target_metrics(FileO):
    for line in FileO:
        if line.startswith('BAIT_SET'):
            break
    reader=csv.DictReader(FileO,delimiter='\t',fieldnames=line.rstrip().split('\t'))
    METRICS=reader.__next__()
    return METRICS

def disambiguate_metrics(FileO):
    #DEPRECIATED
    METRICS={}
    reader=csv.DictReader(FileO,delimiter='\t')
    LINE=reader.__next__()
    #values=dict(zip(file.__next__().rstrip().split('\t'),file.__next__().rstrip().split('\t')))
    total=0.0
    in_keys=['unique species A pairs','unique species B pairs','ambiguous pairs']
    for k in in_keys:
        total+=float(LINE.get(k,0.0))
    F=lambda x: x/total
    out_keys=['PCT_HUMAN','PCT_MOUSE','PCT_AMBIGUOUS']
    for k,v in zip(out_keys,in_keys):
        try:
            METRICS[k]=F(float(LINE.get(v,0.0)))
        except (IndexError,ZeroDivisionError) as e:
            print(f'{e}: in {FileO.name}')
    #print(METRICS)
    return METRICS

def refstats_metrics(FileO):
    reader=csv.DictReader(FileO,delimiter='\t')
    x=next(reader)
    y=next(reader)
    return {'PCT_HUMAN':x['%unambiguousReads'],'PCT_MOUSE':y['%unambiguousReads'],'HUMAN_READS':x['unambiguousReads'],'MOUSE_READS':y['unambiguousReads']}

def run_sample(sample,targets,ref='GRCh37',PDX=False):
    print(sample)
    OUT_METRICS=defaultdict(str)
    #Duplicates
    try:
        dup_file=os.path.abspath(f'metrics/{ref}/{sample}/flagstat.metrics')
        with open(dup_file,'r') as FILE:
            OUT_METRICS.update(flagstat_metrics(FILE))
    except FileNotFoundError as e:
        print(f'{e}: {dup_file}')
    
    #Alignment
    try:
        aln_file=os.path.abspath(f'metrics/{ref}/{sample}/alignment_summary.metrics')
        with open(aln_file,'r') as FILE:
            OUT_METRICS.update(alignment_metrics(FILE))
    except (FileNotFoundError,StopIteration) as e:
        print(f'{e}: {aln_file}')
    
    #Target metrics
    try:
        target_file=os.path.abspath(f'metrics/{targets}/{sample}/target.metrics')
        with open(target_file,'r') as FILE:
            OUT_METRICS.update(target_metrics(FILE))
            OUT_METRICS['PCT_TARGET_BASES_lt_20X']=str(1.0-float(OUT_METRICS['PCT_TARGET_BASES_20X']))
    except FileNotFoundError as e:
        print(f'{e}: {target_file}')
    
    if PDX:
        #clean this up
        print(os.listdir(f"metrics/{sample}/"))
        try:
            refstats_files=glob.glob(f'metrics/{sample}/*.refstats')
            assert len(refstats_files)>0
            x=[]
            y=[]
            for f in refstats_files:
                with open(f,'r') as FILE:
                    m=refstats_metrics(FILE)
                    x.append(m['PCT_HUMAN'])
                    y.append(m['PCT_MOUSE'])
            OUT_METRICS['PCT_HUMAN']=str(round(sum(list(map(float,x)))/len(x)*0.01,3))
            OUT_METRICS['PCT_MOUSE']=str(round(sum(list(map(float,y)))/len(y)*0.01,3))
        except FileNotFoundError as e:
            print(f'{e}: What. {sample}')
    OUT_METRICS['SAMPLE']=sample
    return OUT_METRICS
    #DEPRECIATED
    #Now using bbsplit
    #Disambig
#    if PDX:
#        try:
#            disambres_file=os.path.abspath(f'bam_input/work/{sample}/disambres/mapped_summary.txt')
#            with open(disambres_file,'r') as FILE:
#                OUT_METRICS.update(disambiguate_metrics(FILE))
#        except FileNotFoundError as e:
#            print(f'{e}: {disambres_file}')





def format_outrow(OUT_METRICS):
    FINAL_METRICS={}
    for k,v in OUT_METRICS.items():
        if v=='':
            FINAL_METRICS[k]='.'
        elif 'PCT' in k:
            try:
                FINAL_METRICS[k]=str(round(float(v),3))
            except ValueError:
                FINAL_METRICS[k]=v
        else:
            FINAL_METRICS[k]=str(v)
    return FINAL_METRICS

def write_summary(SAMPLES,output_fp,targets,ref,PDX):
    if PDX==True:
        header=header='SAMPLE BAIT_SET TOTAL_READS PCT_HUMAN PCT_MOUSE PCT_PF_READS PF_HQ_ALIGNED_READS PCT_READS_ALIGNED_IN_PAIRS PCT_DUP PCT_SELECTED_BASES MEAN_TARGET_COVERAGE MEDIAN_TARGET_COVERAGE MAX_TARGET_COVERAGE PCT_USABLE_BASES_ON_TARGET ZERO_CVG_TARGETS_PCT PCT_TARGET_BASES_20X PCT_TARGET_BASES_100X'.split(' ')
    else:
        header='SAMPLE BAIT_SET TOTAL_READS PCT_PF_READS PF_HQ_ALIGNED_READS PCT_READS_ALIGNED_IN_PAIRS PCT_DUP PCT_SELECTED_BASES MEAN_TARGET_COVERAGE MEDIAN_TARGET_COVERAGE MAX_TARGET_COVERAGE PCT_USABLE_BASES_ON_TARGET ZERO_CVG_TARGETS_PCT PCT_TARGET_BASES_20X PCT_TARGET_BASES_100X'.split(' ')
    with open(output_fp,'w') as OFILE:
        writer=csv.DictWriter(OFILE,delimiter=',',fieldnames=header,restval='.',extrasaction='ignore')
        writer.writeheader()
        for sample in SAMPLES:
            SAMPLE_METRICS=format_outrow(run_sample(sample,targets,ref,PDX))
            writer.writerow(SAMPLE_METRICS)

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-I','--input_fp',default='samples.list',help='Sample list')
    p.add_argument('-O','--output_fp',default='metrics_summary.csv',help='Output metrics.csv')
    p.add_argument('-L','--targets',required=True,help='targets key')
    p.add_argument('-R','--ref',default='GRCh37',help='reference key')
    p.add_argument('--PDX',choices=['True','False'],default='False',help='PDX disambiguation from Mouse')
    return p.parse_args()

def main(argv=None):
    args=get_args()
    input=args.input_fp
    output=args.output_fp
    ref=args.ref
    targets=args.targets
    PDX=(args.PDX=='True')
    print(f'-I {input} -O {output} -L {targets} -R {ref} --PDX {PDX}')
    with open(input,'r') as file:
        SAMPLES=file.read().splitlines()
    write_summary(SAMPLES,output,targets,ref,PDX)


if __name__=='__main__':
    main()