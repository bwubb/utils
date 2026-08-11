import os
import csv
import yaml
import argparse

def target_metrics(FileO):
    for line in FileO:
        if line.startswith('BAIT_SET'):
            break
    reader=csv.DictReader(FileO,delimiter='\t',fieldnames=line.rstrip().split('\t'))
    return reader.__next__()

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

def load_pairs(fastq_config_fp,sample_list_fp):
    with open(fastq_config_fp) as file:
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
            if r1 is None or r2 is None:
                continue
            paired.append((base_stem,r1,r2))
        if not paired:
            continue
        new_val={}
        for base_stem,r1,r2 in paired:
            parts=base_stem.split('_')
            run,lane,index=(parts[-3],parts[-2],parts[-1]) if len(parts)>=3 else ('run','0','0')
            new_val[base_stem]={'run':run,'lane':lane,'index':index,'files':[r1,r2]}
        FILES[sample]=new_val
    with open(sample_list_fp) as file:
        SAMPLES=[s for s in file.read().splitlines() if s in FILES]
    pairs=[]
    for sample in SAMPLES:
        for base_stem in sorted(FILES[sample].keys()):
            entry=FILES[sample][base_stem]
            pairs.append((sample,entry['run'],entry['lane'],entry['index']))
    return pairs

def markdup_metrics(FileO):
    stats={}
    for line in FileO:
        if ':' not in line:
            continue
        k,v=line.split(':',1)
        stats[k.strip()]=v.strip()
    read=int(stats['READ'])
    dup=int(stats['DUPLICATE PRIMARY TOTAL'])
    return {'PCT_DUP':str(float(dup/read))}

def run_pair(sample,run,lane,index,ref):
    base=os.path.abspath(f"bam_input/work/{sample}/{ref}/{run}/{lane}/{index}")
    OUT={'SAMPLE':sample,'RUN':run,'LANE':lane,'INDEX':index}
    stats_fp=f"{base}/5.stats.txt"
    hs_fp=f"{base}/5.hs.metrics"
    try:
        with open(stats_fp,'r') as FILE:
            OUT.update(markdup_metrics(FILE))
    except (FileNotFoundError,KeyError,ValueError,ZeroDivisionError) as e:
        print(f'{e}: {stats_fp}')
    try:
        with open(hs_fp,'r') as FILE:
            hs=target_metrics(FILE)
            OUT['TOTAL_READS']=hs.get('TOTAL_READS','')
            OUT['MEAN_TARGET_COVERAGE']=hs.get('MEAN_TARGET_COVERAGE','')
            OUT['PCT_USABLE_BASES_ON_TARGET']=hs.get('PCT_USABLE_BASES_ON_TARGET','')
    except (FileNotFoundError,StopIteration) as e:
        print(f'{e}: {hs_fp}')
    return OUT

def write_summary(pairs,ref,output_fp):
    header='SAMPLE RUN LANE INDEX PCT_DUP TOTAL_READS MEAN_TARGET_COVERAGE PCT_USABLE_BASES_ON_TARGET'.split(' ')
    with open(output_fp,'w') as OFILE:
        writer=csv.DictWriter(OFILE,delimiter=',',fieldnames=header,restval='.',extrasaction='ignore')
        writer.writeheader()
        for sample,run,lane,index in pairs:
            writer.writerow(format_outrow(run_pair(sample,run,lane,index,ref)))

def get_args(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('-I','--input_fp',default='sample.list',help='Sample list')
    p.add_argument('-F','--fastq_config',default='fastq.yml',help='FASTQ config')
    p.add_argument('-O','--output_fp',default='pair_metrics_summary.csv',help='Output CSV')
    p.add_argument('-R','--ref',default='GRCh38',help='Reference key')
    return p.parse_args(argv)

def main(argv=None):
    args=get_args(argv)
    pairs=load_pairs(args.fastq_config,args.input_fp)
    write_summary(pairs,args.ref,args.output_fp)

if __name__=='__main__':
    main()
