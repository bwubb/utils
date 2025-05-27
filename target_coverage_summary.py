import os
import csv
import argparse
from collections import defaultdict

def get_args(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('-I','--input_fp',required=True,help='Sample list')
    p.add_argument('-O','--output_fp',required=True,help='Output coverage summary CSV')
    p.add_argument('-L','--targets',required=True,help='targets key')
    return p.parse_args(argv)

def main(argv=None):
    args=get_args(argv)
    with open(args.input_fp,'r') as file:
        SAMPLES=file.read().splitlines()
    
    target_coverage=defaultdict(lambda: defaultdict(float))
    key_order=[]
    
    for sample in SAMPLES:
        sample_file=os.path.abspath(f"metrics/{args.targets}/{sample}/target_coverage.metrics")
        with open(sample_file,'r') as cov_file:
            reader=csv.DictReader(cov_file,delimiter='\t')
            for row in reader:
                key=tuple([row[x] for x in ['chrom','start','end','length','name']])
                if key not in key_order:
                    key_order.append(key)
                target_coverage[key][sample]=float(row['mean_coverage'])
    
    with open(args.output_fp,'w') as outfile:
        writer=csv.writer(outfile,delimiter=',')
        writer.writerow(['Chr','Start','End','Length','Name']+SAMPLES)
        for key in key_order:
            row=list(key)+[target_coverage[key][sample] for sample in SAMPLES]
            writer.writerow(row)

if __name__=='__main__':
    main() 