#!/usr/bin/env python

import argparse
import datetime
import glob
import csv
import os
import shutil
import logging
from Bio.Seq import Seq

def rev_comp(seq):
    return str(Seq(seq).reverse_complement())

def check_index_variants(index1,index2):
    variants=[]
    if index1 and index2:
        variants.append(f'{index1}-{index2}')
        variants.append(f'{index2}-{index1}')
        variants.append(f'{rev_comp(index1)}-{index2}')
        variants.append(f'{index1}-{rev_comp(index2)}')
        variants.append(f'{rev_comp(index1)}-{rev_comp(index2)}')
    return variants

def find_FGC_fastqs(dir,run,lane,index_variants,sample):
    for index in index_variants:
        fqs=sorted(glob.glob(f'{dir}/{run}_s_{lane}_[1-2]_{index}.fastq.gz'))
        if len(fqs)==2:
            return fqs,index
    return [],None

def find_Nextseq_fastqs(dir,run,lane,index_variants,sample):
    for index in index_variants:
        fqs=sorted(glob.glob(f'{dir}/{sample}_S[0-9]*_L00{lane}_R[1-2]_001.fastq.gz'))
        if len(fqs)==2:
            return fqs,index
    return [],None

def validate_input_file(infile):
    required_cols=['sample','run','lane','index','lib']
    with open(infile,'r') as f:
        reader=csv.reader(f,delimiter='\t')
        headers=next(reader)
        headers_lower=[h.lower() for h in headers]
        missing=[col for col in required_cols if not any(col in h for h in headers_lower)]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")
        return headers

def main(argv=None):
    logging.basicConfig(level=logging.INFO,format='%(message)s')
    logger=logging.getLogger(__name__)
    
    try:
        headers=validate_input_file(argv.infile)
    except ValueError as e:
        logger.error(f"Input file validation failed: {e}")
        return
    
    m=list(os.path.splitext(os.path.basename(argv.infile)))
    m.insert(-1,'.report')
    
    stats={'Found':0,'Found (Index Variant)':0,'Missing':0,'Error':0,'Already Exists':0,'Unexpected':0}
    
    with open(argv.infile,'r') as infile,open(''.join(m),'w') as report:
        reader=csv.reader(infile,delimiter='\t')
        next(reader)
        report_writer=csv.writer(report,delimiter='\t')
        
        # Write header
        report_writer.writerow(['Fastq Processing Report'])
        report_writer.writerow(['Date:',datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
        report_writer.writerow(['Input:',argv.infile])
        report_writer.writerow(['Dir:',argv.dir])
        report_writer.writerow(['Mode:',argv.action])
        report_writer.writerow([])
        
        report_writer.writerow(['Results'])
        report_writer.writerow(['Sample','Run','Lane','Expected','Found','Status','Note'])
        report_writer.writerow([])  # Empty line instead of dashes
        
        sample=False
        run=False
        lane=False
        lib=False
        indices=[]
        for field in headers:
            if 'sample' in field.lower():
                sample=headers.index(field)
            elif 'run' in field.lower():
                run=headers.index(field)
            elif 'lane' in field.lower():
                lane=headers.index(field)
            elif any([x in field.lower() for x in ['index','barcode']]):
                indices.append(headers.index(field))
            elif 'lib' in field.lower():
                lib=headers.index(field)
        assert not all(v is False for v in [run,lane,sample,len(indices)>0])
        
        if argv.source=='fgc':
            f=find_FGC_fastqs
        elif argv.source=='nextseq':
            f=find_Nextseq_fastqs
        else:
            raise
        
        found_files=set()
        for row in reader:
            expected_index="-".join([row[x] for x in indices])
            index_variants=check_index_variants(*[row[x] for x in indices])
            fqs,found_index=f(argv.dir,row[run],row[lane],index_variants,row[sample])
            
            new_R1_name=f'{row[sample]}_{row[lib]}_{row[run]}_{row[lane]}_{expected_index}_R1.fastq.gz'
            new_R2_name=f'{row[sample]}_{row[lib]}_{row[run]}_{row[lane]}_{expected_index}_R2.fastq.gz'
            
            if len(fqs)==2:
                found_files.add(fqs[0])
                found_files.add(fqs[1])
                try:
                    logger.info(f"Found: {fqs[0]} -> {new_R1_name}")
                    logger.info(f"Found: {fqs[1]} -> {new_R2_name}")
                    
                    if argv.action=='rename':
                        os.rename(fqs[0],f'{argv.dir}/{new_R1_name}')
                        os.rename(fqs[1],f'{argv.dir}/{new_R2_name}')
                    elif argv.action=='copy':
                        shutil.copy2(fqs[0],f'{argv.dir}/{new_R1_name}')
                        shutil.copy2(fqs[1],f'{argv.dir}/{new_R2_name}')
                    
                    status='Found'
                    notes=''
                    if found_index!=expected_index:
                        status='Found (Index Variant)'
                        notes=f'Expected {expected_index}, found {found_index}'
                    
                    stats[status]+=1
                    report_writer.writerow([row[sample],row[run],row[lane],expected_index,found_index,status,notes])
                except (OSError,IOError) as e:
                    logger.error(f"Error processing files: {e}")
                    stats['Error']+=1
                    report_writer.writerow([row[sample],row[run],row[lane],expected_index,found_index,'Error',str(e)])
            elif os.path.isfile(f'{argv.dir}/{new_R1_name}'):
                logger.info(f'File already exists: {new_R1_name}')
                stats['Already Exists']+=1
                report_writer.writerow([row[sample],row[run],row[lane],expected_index,expected_index,'Already Exists',''])
            else:
                logger.warning(f'No files found for {row[sample]}')
                stats['Missing']+=1
                report_writer.writerow([row[sample],row[run],row[lane],expected_index,'','Missing',''])
        
        all_fastqs=set(glob.glob(f'{argv.dir}/*.fastq.gz'))
        unexpected=all_fastqs-found_files
        if unexpected:
            report_writer.writerow(['','','','','','Unexpected Files Found',''])
            for fq in sorted(unexpected):
                report_writer.writerow(['','','','',os.path.basename(fq),'Unexpected',''])

    # Write summary after file is closed
    with open(''.join(m),'a') as report:
        report_writer=csv.writer(report,delimiter='\t')
        report_writer.writerow([])
        report_writer.writerow(['Summary'])
        for status,count in stats.items():
            if count>0:
                report_writer.writerow([f'{status}: {count}'])
        
        if unexpected:
            report_writer.writerow([])
            report_writer.writerow(['Extra Files'])
            report_writer.writerow([])  # Empty line instead of dashes
            for fq in sorted(unexpected):
                report_writer.writerow([os.path.basename(fq)])

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('-i','--infile',help='Sequence Project submission file')
    p.add_argument('-D','--dir',default='./FASTQ',help='Fastq directory')
    p.add_argument('--source',choices=['fgc','nextseq'],default='fgc',help='Was this from the FGC or NextSeq?')
    p.add_argument('--action',choices=['copy','rename','dryrun'],default='dryrun',help='What action to take.')
    argv=p.parse_args()
    logger=logging.getLogger(__name__)
    logger.info("Arguments selected:")
    for a,b in vars(argv).items():
        logger.info(f' {a}: {b}')
    main(argv)
