
import yaml
import glob
import re
import os
import datetime
import argparse
from collections import defaultdict
from yaml.representer import Representer
yaml.add_representer(defaultdict, Representer.represent_dict)

def subset_R2_fastqs(FASTQS):
    r2=[]
    for file in FASTQS:
        if re.search('.*_R2\.fastq\.gz',os.path.basename(file)):
            r2.append(os.path.basename(file))
    return r2

def make_info(r2):
    sample_data=defaultdict(dict)
    for fq2 in r2:#need os basename
        basename=fq2.rstrip('.fastq.gz')
        parts=basename.split('_')
        #assert len
        sample_name=parts[0]
        LB=parts[1]
        ID='-'.join(parts[2:4])
        PU='-'.join(parts[2:5])
        fq1=fq2.replace('_R2.fastq.gz','_R1.fastq.gz')
        #if files exist else ~
        files=[fq1,fq2]
        sample_data[sample_name][ID]={'LB':LB,'PU':PU,'files':files}
    return sample_data

def write_output(outdict,filename):
    with open(f'fastq.yaml','w') as outfile:
        yaml.dump(outdict,outfile)
    with open(f'sample.list','w') as file:
        for i in sorted(outdict.keys()):
            file.write(f'{i}\n')

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('--dir',default='FASTQ',help='FASTQ dir')
    p.add_argument('--date',action='store_true',default=False,help='add date to files')
    args=p.parse_args()
    
    all_fastqs=glob.glob(f'{args.dir}/*.fastq.gz')
    r2=subset_R2_fastqs(all_fastqs)
    sample_data=make_info(r2)
    
    validate=defaultdict(set)
    for sample_key,sample_vals in sample_data.items():
        for runid_key,runid_vals in sample_vals.items():
            validate[sample_key].add(runid_vals['LB'])
    
    with open(f'{datetime.date.today().strftime("%Y%m%d")}.sample.targets.yaml','w') as file:
        #Custom YAML write
        #Not sorted
        for sample_key,LB in validate.items():
            if len(LB)!=1:
                print('WARNING! Multiple LB values for {sample_key} : {",".join(LB)}')
                #Write possible fix, to rename files, and sample names with *-2,*-3,etc
            else:
                file.write(f'{sample_key}: {"".join(LB)}\n')
    
    write_output(sample_data,f'{datetime.date.today().strftime("%Y%m%d")}')
    
    outfile_dict=defaultdict(lambda: defaultdict(dict))
    for Sample,v in sample_data.items():
        for RunLane,v1 in v.items():
            if RunLane not in outfile_dict[v1['LB']][Sample].keys():
                outfile_dict[v1['LB']][Sample][RunLane]=v1
    for key in outfile_dict.keys():
        write_output(outfile_dict[key],f'{datetime.date.today().strftime("%Y%m%d")}.{key}')


if __name__=='__main__':
    main()

#full project.yaml with sample id, and lib and bam and ref version of the bam, fasta, etc.
# resource yaml
