
import yaml
import glob
import re
import os
import datetime
from collections import defaultdict
from yaml.representer import Representer
yaml.add_representer(defaultdict, Representer.represent_dict)

def subset_R2_fastqs(FASTQS):
    r2=[]
    for file in FASTQS:
        #print(file)
        if re.search('.*_R2\.fastq\.gz',os.path.basename(file)):
            r2.append(os.path.basename(file))
        else:
            pass
    print(len(r2))
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
    
    
def main(argv=None):
    all_fastqs=glob.glob('FASTQ/*.fastq.gz')
    #print(all_fastqs)
    r2=subset_R2_fastqs(all_fastqs)
    sample_data=make_info(r2)


    with open(f'fastqs.{datetime.date.today().strftime("%Y%m%d")}.yaml','w') as outfile:
        yaml.dump(sample_data,outfile)

    with open(f'samples.{datetime.date.today().strftime("%Y%m%d")}.list','w') as file:
        for sample in sample_data.keys():
            file.write(f'{sample}\n')

#    for key in sample_data.keys():
#        print(key)
#        with open(f'bam_input/final/{key}/GRCh37/fastqs.yaml','w') as outfile:
#            yaml.dump(sample_data[key],outfile)

if __name__=='__main__':
    main()
