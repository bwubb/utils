#!/usr/bin/env python

import argparse
import datetime
import glob
import csv
import os

def find_FGC_fastqs(dir,run,lane,index):
    fqs=sorted(glob.glob(f'{dir}/{run}_s_{lane}_[1-2]_{index}.fastq.gz'))
    if len(fqs)!=2:
        return None
    else:
        return fqs


def main(argv=None):
    m=list(os.path.splitext(os.path.basename(argv.infile)))
    m.insert(-1,'.missing')
    with open(argv.infile,'r') as infile,open(''.join(m),'w') as missing:
        reader=csv.reader(infile,delimiter='\t')
        fields=next(reader)
        writer=csv.DictWriter(missing,delimiter='\t',fieldnames=fields)
        writer.writeheader()
        sample=False
        run=False
        lane=False
        lib=argv.lib
        indices=[]
        for field in fields:
            print(field)
            if 'sample' in field.lower():
                sample=fields.index(field)
            elif 'run' in field.lower():
                run=fields.index(field)
            elif 'lane' in field.lower():
                lane=fields.index(field)
            elif any([x in field.lower() for x in ['index','barcode']]):
                indices.append(fields.index(field))
            

        #This is terrible.
        assert not all(v is False for v in [run,lane,sample,len(indices)>0])#I have found the columns
        
        #CLEAN THIS ALLL UP YOU AMATEUR
        for row in reader:
            index="-".join([row[x] for x in indices])
            fqs=find_FGC_fastqs(argv.dir,row[run],row[lane],index)
            new_R1_name=f'{row[sample]}_{lib}_{row[run]}_{row[lane]}_{index}_R1.fastq.gz'
            new_R2_name=f'{row[sample]}_{lib}_{row[run]}_{row[lane]}_{index}_R2.fastq.gz'
            if fqs:
                print(fqs[0],f'> {argv.dir}/{new_R1_name}')
                print(fqs[1],f'> {argv.dir}/{new_R2_name}')
                if argv.action=='rename':
                    os.rename(fqs[0],f'{argv.dir}/{new_R1_name}')
                    os.rename(fqs[1],f'{argv.dir}/{new_R2_name}')
        #actions

if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('-i','--infile',help='Sequence Project submission file')
    p.add_argument('-D','--dir',default='./FASTQ',help='Fastq directory')
    p.add_argument('-L','--lib',help='Library Targets Key')
    p.add_argument('--action',choices=['copy','rename','dryrun'],default='dryrun',help='What action to take.')
    argv=p.parse_args()
    print("Arguments selected:")
    for a, b in vars(argv).items():
        print(f' {a}: {b}')
    main(argv)