import os
import re
from pathlib import Path
from .config import REFERENCE_PATHS

class IntervalManager:
    def __init__(self):
        self.chr_pattern=re.compile(r'^(chr)?(.+)$')
    
    def _handle_chr_prefix(self,chrom,add_prefix):
        """Add or remove chr prefix from chromosome name"""
        match=self.chr_pattern.match(chrom)
        if not match:
            return chrom
        
        has_prefix,rest=match.groups()
        if add_prefix and not has_prefix:
            return f'chr{rest}'
        elif not add_prefix and has_prefix:
            return rest
        return chrom
    
    def _load_chromosome_lengths(self,dict_file):
        """Load chromosome lengths from reference dictionary file"""
        chrom_lengths={}
        with open(dict_file,'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    parts=line.strip().split('\t')
                    chrom=None
                    length=None
                    for part in parts[1:]:
                        if part.startswith('SN:'):
                            chrom=part[3:]
                        elif part.startswith('LN:'):
                            length=int(part[3:])
                    if chrom and length:
                        chrom_lengths[chrom]=length
        return chrom_lengths
    
    def validate_bed(self,args):
        """Validate BED file coordinates against reference dictionary"""
        chrom_lengths=self._load_chromosome_lengths(args.dict)
        
        with open(args.input,'r') as bed,open(args.output,'w') as out:
            for line_num,line in enumerate(bed,1):
                if line.startswith('#'):
                    out.write(line)
                    continue
                
                fields=line.strip().split('\t')
                if len(fields)<3:
                    continue
                
                chrom=fields[0]
                start=int(fields[1])
                end=int(fields[2])
                
                # Check if chromosome exists in dict
                if chrom not in chrom_lengths:
                    print(f"Warning: Chromosome {chrom} not found in dict file, skipping line {line_num}")
                    continue
                
                # Check if coordinates are within bounds
                if start<0 or end>chrom_lengths[chrom] or start>=end:
                    print(f"Warning: Invalid coordinates on {chrom}:{start}-{end} (length: {chrom_lengths[chrom]}), skipping line {line_num}")
                    continue
                
                out.write(line)
    
    def create_picard_intervals(self,args):
        """Create Picard intervals from bait and target files"""
        chrom_lengths=self._load_chromosome_lengths(args.dict)
        
        # Get output filenames
        bait_basename=os.path.splitext(args.bait)[0]
        target_basename=os.path.splitext(args.target)[0]
        bait_output=f"{bait_basename}.picard_bait.intervals"
        target_output=f"{target_basename}.picard_target.intervals"
        
        # Validate chromosomes in bait file
        bait_chroms=set()
        with open(args.bait,'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields=line.strip().split('\t')
                if len(fields)>=3:
                    bait_chroms.add(fields[0])
        
        # Validate chromosomes in target file
        target_chroms=set()
        with open(args.target,'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields=line.strip().split('\t')
                if len(fields)>=3:
                    target_chroms.add(fields[0])
        
        # Check for missing chromosomes in dict
        dict_chroms=set(chrom_lengths.keys())
        missing_bait=bait_chroms-dict_chroms
        missing_target=target_chroms-dict_chroms
        
        if missing_bait:
            print(f"Warning: Bait chromosomes not in dict: {missing_bait}")
        if missing_target:
            print(f"Warning: Target chromosomes not in dict: {missing_target}")
        
        # Process bait file
        with open(bait_output,'w') as out:
            with open(args.bait,'r') as bait:
                for line in bait:
                    if line.startswith('#'):
                        continue
                    fields=line.strip().split('\t')
                    if len(fields)<3:
                        continue
                    
                    chrom=fields[0]
                    start=int(fields[1])
                    end=int(fields[2])
                    
                    # Validate coordinates
                    if chrom not in chrom_lengths or start<0 or end>chrom_lengths[chrom] or start>=end:
                        continue
                    
                    # Use strand from field 5 if available, otherwise '+'
                    strand=fields[5] if len(fields)>5 else '+'
                    # Use name from field 3 if available, otherwise generate one
                    name=fields[3] if len(fields)>3 else f"{chrom}_{start}_{end}"
                    
                    out.write(f"{chrom}\t{start}\t{end}\t{strand}\t{name}\n")
        
        # Process target file
        with open(target_output,'w') as out:
            with open(args.target,'r') as target:
                for line in target:
                    if line.startswith('#'):
                        continue
                    fields=line.strip().split('\t')
                    if len(fields)<3:
                        continue
                    
                    chrom=fields[0]
                    start=int(fields[1])
                    end=int(fields[2])
                    
                    # Validate coordinates
                    if chrom not in chrom_lengths or start<0 or end>chrom_lengths[chrom] or start>=end:
                        continue
                    
                    # Use strand from field 5 if available, otherwise '+'
                    strand=fields[5] if len(fields)>5 else '+'
                    # Use name from field 3 if available, otherwise generate one
                    name=fields[3] if len(fields)>3 else f"{chrom}_{start}_{end}"
                    
                    out.write(f"{chrom}\t{start}\t{end}\t{strand}\t{name}\n")
    
    def bed_to_intervals(self,bed_file,output_file,add_chr_prefix=False):
        """Convert BED to GATK intervals format"""
        with open(bed_file,'r') as bed,open(output_file,'w') as out:
            for line in bed:
                if line.startswith('#'):
                    continue
                fields=line.strip().split('\t')
                if len(fields)<3:
                    continue
                
                chrom=self._handle_chr_prefix(fields[0],add_chr_prefix)
                start=int(fields[1])+1  # BED is 0-based, intervals are 1-based
                end=int(fields[2])
                out.write(f'{chrom}:{start}-{end}\n')
    
    def intervals_to_bed(self,intervals_file,output_file,add_chr_prefix=False):
        """Convert GATK intervals to BED format"""
        with open(intervals_file,'r') as intervals,open(output_file,'w') as out:
            for line in intervals:
                if line.startswith('#'):
                    continue
                chrom,pos=line.strip().split(':')
                chrom=self._handle_chr_prefix(chrom,add_chr_prefix)
                start,end=map(int,pos.split('-'))
                out.write(f'{chrom}\t{start-1}\t{end}\n')  # Convert to 0-based
    
    def bed_to_bed(self,bed_file,output_file,add_chr_prefix=False):
        """Convert BED to BED format with chr prefix handling"""
        with open(bed_file,'r') as bed,open(output_file,'w') as out:
            for line in bed:
                if line.startswith('#'):
                    out.write(line)
                    continue
                fields=line.strip().split('\t')
                if len(fields)<3:
                    continue
                chrom=self._handle_chr_prefix(fields[0],add_chr_prefix)
                # Preserve all columns from input
                fields[0]=chrom
                out.write('\t'.join(fields)+'\n')
    
    def intervals_to_intervals(self,intervals_file,output_file,add_chr_prefix=False):
        """Convert intervals to intervals format with chr prefix handling"""
        with open(intervals_file,'r') as intervals,open(output_file,'w') as out:
            for line in intervals:
                if line.startswith('#'):
                    out.write(line)
                    continue
                chrom,pos=line.strip().split(':')
                chrom=self._handle_chr_prefix(chrom,add_chr_prefix)
                out.write(f'{chrom}:{pos}\n')
    
    def make_picard_intervals(self,bed_file,output_file,reference,add_chr_prefix=False):
        """Create Picard-style intervals file"""
        # Get the reference dictionary file
        dict_file=REFERENCE_PATHS[reference]['dict']
        if not os.path.exists(dict_file):
            raise FileNotFoundError(f"Reference dictionary file not found: {dict_file}")
        
        # First read the dict file to get chromosome lengths
        chrom_lengths={}
        with open(dict_file,'r') as f:
            for line in f:
                if line.startswith('@SQ'):
                    fields=dict(f.split(':') for f in line.strip().split('\t')[1:])
                    chrom=fields['SN']
                    length=int(fields['LN'])
                    chrom_lengths[chrom]=length
        
        # Now write the intervals
        with open(bed_file,'r') as bed,open(output_file,'w') as out:
            # Write the dict header
            for chrom,length in chrom_lengths.items():
                out.write(f'@SQ\tSN:{chrom}\tLN:{length}\n')
            
            # Write the intervals
            for line in bed:
                if line.startswith('#'):
                    continue
                fields=line.strip().split('\t')
                if len(fields)<3:
                    continue
                
                chrom=self._handle_chr_prefix(fields[0],add_chr_prefix)
                start=int(fields[1])+1  # BED is 0-based
                end=int(fields[2])
                
                out.write(f'{chrom}\t{start}\t{end}\n')
    
    def process_intervals(self,args):
        """Main entry point for interval processing"""
        add_chr_prefix=args.reference=='hg38'
        
        # Check if input is BED or intervals format
        is_bed=True
        with open(args.input,'r') as f:
            first_line=f.readline().strip()
            if ':' in first_line and '-' in first_line:
                is_bed=False
        
        if args.format=='bed':
            if is_bed:
                self.bed_to_bed(args.input,args.output,add_chr_prefix)
            else:
                self.intervals_to_bed(args.input,args.output,add_chr_prefix)
        elif args.format=='intervals':
            if is_bed:
                self.bed_to_intervals(args.input,args.output,add_chr_prefix)
            else:
                self.intervals_to_intervals(args.input,args.output,add_chr_prefix) 