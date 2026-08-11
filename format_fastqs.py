
import yaml
import glob
import re
import os
import datetime
import argparse
from collections import defaultdict
from yaml.representer import Representer
yaml.add_representer(defaultdict, Representer.represent_dict)

# fastq.yml for alignment.smk: flat per sample
#   SAMPLE:
#     files:
#       - ..._R1.fastq.gz
#       - ..._R2.fastq.gz
# (basenames; alignment.smk prepends FASTQ/ unless path is absolute)

def subset_R2_fastqs(FASTQS):
    r2=[]
    for file in FASTQS:
        if re.search('.*_R2\.fastq\.gz',os.path.basename(file)):
            r2.append(os.path.basename(file))
    return r2

def make_info(r2):
    sample_data=defaultdict(dict)
    for fq2 in r2:
        basename=fq2.rstrip('.fastq.gz')
        parts=basename.split('_')
        sample_name=parts[0]
        LB=parts[1]
        run,lane,index=(parts[-3],parts[-2],parts[-1]) if len(parts)>=3 else ('run','0','0')
        fq1=fq2.replace('_R2.fastq.gz','_R1.fastq.gz')
        files=[fq1,fq2]
        sample_data[sample_name][basename]={'LB':LB,'run':run,'lane':lane,'index':index,'files':files}
    return sample_data

def flat_fastq_config(sample_data):
    """Collapse nested run entries to sample -> files (alignment.smk / project-utils shape)."""
    flat={}
    for sample,runs in sample_data.items():
        files=[]
        for entry in runs.values():
            files.extend(entry['files'])
        flat[sample]={'files':sorted(files)}
    return flat

def write_output(flat_config):
    with open('fastq.yml','w') as outfile:
        yaml.dump(flat_config,outfile,default_flow_style=False)
    with open('sample.list','w') as file:
        for i in sorted(flat_config.keys()):
            file.write(f'{i}\n')
    print('Created: fastq.yml, sample.list')

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('--dir',default='FASTQ',help='FASTQ dir')
    p.add_argument('--date',action='store_true',default=False,help='add date to files')
    args=p.parse_args()

    all_fastqs=glob.glob(f'{args.dir}/*.fastq.gz')
    print(f"{len(all_fastqs)} found.")
    r2=subset_R2_fastqs(all_fastqs)
    sample_data=make_info(r2)

    validate=defaultdict(set)
    for sample_key,sample_vals in sample_data.items():
        for runid_key,runid_vals in sample_vals.items():
            validate[sample_key].add(runid_vals['LB'])

    date_prefix=datetime.date.today().strftime('%Y%m%d')
    with open(f'{date_prefix}.sample.targets.yml','w') as file:
        for sample_key,LB in validate.items():
            if len(LB)!=1:
                print(f'WARNING! Multiple LB values for {sample_key}: {",".join(LB)}')
            else:
                file.write(f'{sample_key}: {"".join(LB)}\n')

    write_output(flat_fastq_config(sample_data))

if __name__=='__main__':
    main()

#full project.yaml with sample id, and lib and bam and ref version of the bam, fasta, etc.
# resource yaml
