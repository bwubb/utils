import os
import yaml
import glob
import re
import csv
import shutil
import logging
import datetime
from pathlib import Path
from collections import defaultdict
from Bio.Seq import Seq

class SampleManager:
    def __init__(self):
        self.logger=logging.getLogger(__name__)
    
    def rev_comp(self,seq):
        """Get reverse complement of sequence"""
        return str(Seq(seq).reverse_complement())
    
    def check_index_variants(self,index1,index2):
        """Generate all possible index combinations"""
        variants=[]
        if index1 and index2:
            variants.append(f'{index1}-{index2}')
            variants.append(f'{index2}-{index1}')
            variants.append(f'{self.rev_comp(index1)}-{index2}')
            variants.append(f'{index1}-{self.rev_comp(index2)}')
            variants.append(f'{self.rev_comp(index1)}-{self.rev_comp(index2)}')
        return variants
    
    def find_FGC_fastqs(self,dir,run,lane,index_variants,sample):
        """Find FASTQ files from FGC"""
        for index in index_variants:
            fqs=sorted(glob.glob(f'{dir}/{run}_s_{lane}_[1-2]_{index}.fastq.gz'))
            if len(fqs)==2:
                return fqs,index
        return [],None
    
    def find_Nextseq_fastqs(self,dir,run,lane,index_variants,sample):
        """Find FASTQ files from NextSeq"""
        for index in index_variants:
            fqs=sorted(glob.glob(f'{dir}/{sample}_S[0-9]*_L00{lane}_R[1-2]_001.fastq.gz'))
            if len(fqs)==2:
                return fqs,index
        return [],None
    
    def validate_input_file(self,infile):
        """Validate submission file has required columns"""
        required_cols=['sample','run','lane','index','lib']
        with open(infile,'r') as f:
            reader=csv.reader(f,delimiter='\t')
            headers=next(reader)
            headers_lower=[h.lower() for h in headers]
            missing=[col for col in required_cols if not any(col in h for h in headers_lower)]
            if missing:
                raise ValueError(f"Missing required columns: {', '.join(missing)}")
            return headers
    
    def process_submission_file(self,args):
        """Process sample submission file"""
        try:
            headers=self.validate_input_file(args.input)
        except ValueError as e:
            self.logger.error(f"Input file validation failed: {e}")
            return
        
        # Create report file
        report_file=f"{os.path.splitext(os.path.basename(args.input))[0]}.report.txt"
        stats={'Found':0,'Found (Index Variant)':0,'Missing':0,'Error':0,'Already Exists':0,'Unexpected':0}
        
        with open(args.input,'r') as infile,open(report_file,'w') as report:
            reader=csv.reader(infile,delimiter='\t')
            next(reader)  # Skip header
            report_writer=csv.writer(report,delimiter='\t')
            
            # Write report header
            report_writer.writerow(['Fastq Processing Report'])
            report_writer.writerow(['Date:',datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
            report_writer.writerow(['Input:',args.input])
            report_writer.writerow(['Dir:',args.dir])
            report_writer.writerow(['Mode:',args.action])
            report_writer.writerow([])
            
            report_writer.writerow(['Results'])
            report_writer.writerow(['Sample','Run','Lane','Expected','Found','Status','Note'])
            report_writer.writerow([])
            
            # Find column indices
            sample=run=lane=lib=False
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
            
            # Select finder function
            finder=self.find_FGC_fastqs if args.source=='fgc' else self.find_Nextseq_fastqs
            
            found_files=set()
            for row in reader:
                expected_index="-".join([row[x] for x in indices])
                index_variants=self.check_index_variants(*[row[x] for x in indices])
                fqs,found_index=finder(args.dir,row[run],row[lane],index_variants,row[sample])
                
                new_R1_name=f'{row[sample]}_{row[lib]}_{row[run]}_{row[lane]}_{expected_index}_R1.fastq.gz'
                new_R2_name=f'{row[sample]}_{row[lib]}_{row[run]}_{row[lane]}_{expected_index}_R2.fastq.gz'
                
                if len(fqs)==2:
                    found_files.add(fqs[0])
                    found_files.add(fqs[1])
                    try:
                        self.logger.info(f"Found: {fqs[0]} -> {new_R1_name}")
                        self.logger.info(f"Found: {fqs[1]} -> {new_R2_name}")
                        
                        if args.action=='rename':
                            os.rename(fqs[0],f'{args.dir}/{new_R1_name}')
                            os.rename(fqs[1],f'{args.dir}/{new_R2_name}')
                        elif args.action=='copy':
                            shutil.copy2(fqs[0],f'{args.dir}/{new_R1_name}')
                            shutil.copy2(fqs[1],f'{args.dir}/{new_R2_name}')
                        
                        status='Found'
                        notes=''
                        if found_index!=expected_index:
                            status='Found (Index Variant)'
                            notes=f'Expected {expected_index}, found {found_index}'
                        
                        stats[status]+=1
                        report_writer.writerow([row[sample],row[run],row[lane],expected_index,found_index,status,notes])
                    except (OSError,IOError) as e:
                        self.logger.error(f"Error processing files: {e}")
                        stats['Error']+=1
                        report_writer.writerow([row[sample],row[run],row[lane],expected_index,found_index,'Error',str(e)])
                elif os.path.isfile(f'{args.dir}/{new_R1_name}'):
                    self.logger.info(f'File already exists: {new_R1_name}')
                    stats['Already Exists']+=1
                    report_writer.writerow([row[sample],row[run],row[lane],expected_index,expected_index,'Already Exists',''])
                else:
                    self.logger.warning(f'No files found for {row[sample]}')
                    stats['Missing']+=1
                    report_writer.writerow([row[sample],row[run],row[lane],expected_index,'','Missing',''])
            
            # Check for unexpected files
            all_fastqs=set(glob.glob(f'{args.dir}/*.fastq.gz'))
            unexpected=all_fastqs-found_files
            if unexpected:
                report_writer.writerow(['','','','','','Unexpected Files Found',''])
                for fq in sorted(unexpected):
                    report_writer.writerow(['','','','',os.path.basename(fq),'Unexpected',''])
        
        # Write summary
        with open(report_file,'a') as report:
            report_writer=csv.writer(report,delimiter='\t')
            report_writer.writerow([])
            report_writer.writerow(['Summary'])
            for status,count in stats.items():
                if count>0:
                    report_writer.writerow([f'{status}: {count}'])
            
            if unexpected:
                report_writer.writerow([])
                report_writer.writerow(['Extra Files'])
                report_writer.writerow([])
                for fq in sorted(unexpected):
                    report_writer.writerow([os.path.basename(fq)])
    
    def find_fastq_files(self,input_dir,source='fgc'):
        """Find and organize FASTQ files"""
        fastq_files=defaultdict(list)
        
        if source=='fgc':
            # FGC naming pattern
            pattern='*_R[12]_*.fastq.gz'
        else:
            # NextSeq naming pattern
            pattern='*_S*_R[12]_*.fastq.gz'
        
        for f in glob.glob(os.path.join(input_dir,'**',pattern),recursive=True):
            # Parse sample info from filename
            if source=='fgc':
                sample=os.path.basename(f).split('_R')[0]
            else:
                sample=os.path.basename(f).split('_S')[0]
            fastq_files[sample].append(f)
        
        return fastq_files
    
    def update_sample_list(self,fastq_files,output_dir):
        """Update or create sample list"""
        os.makedirs(output_dir,exist_ok=True)
        
        # Write sample list
        sample_list_file=os.path.join(output_dir,'sample.list')
        with open(sample_list_file,'w') as f:
            for sample in sorted(fastq_files.keys()):
                f.write(f'{sample}\n')
        print(f"Created: {sample_list_file}")
        
        # Write FASTQ config
        fastq_config={'samples':{}}
        for sample,files in fastq_files.items():
            fastq_config['samples'][sample]={
                'files':sorted(files)
            }
        
        fastq_yml_file=os.path.join(output_dir,'fastq.yml')
        with open(fastq_yml_file,'w') as f:
            yaml.dump(fastq_config,f,default_flow_style=False)
        print(f"Created: {fastq_yml_file}")
    
    def find_renamed_fastqs(self,input_dir):
        """Find renamed FASTQ files and organize by sample"""
        fastq_files=defaultdict(list)
        all_fastqs=glob.glob(os.path.join(input_dir,'*.fastq.gz'))
        
        for f in all_fastqs:
            basename=os.path.basename(f)
            # Skip files that start with FGC (original files)
            if not basename.startswith('FGC'):
                # Extract sample name from first part of filename
                sample=basename.split('_')[0]
                fastq_files[sample].append(f)
        
        return fastq_files
    
    def update_samples(self,args):
        """Create sample lists from FASTQ directory"""
        print(f"Scanning directory: {args.dir}")
        fastq_files=self.find_renamed_fastqs(args.dir)
        print(f"Found {len(fastq_files)} samples: {list(fastq_files.keys())}")
        self.update_sample_list(fastq_files,args.dir)
        print(f"Created sample files in: {args.dir}") 