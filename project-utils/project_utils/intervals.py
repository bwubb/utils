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
        elif args.format=='picard':
            if not is_bed:
                # First convert intervals to BED
                temp_bed=args.output+'.temp.bed'
                self.intervals_to_bed(args.input,temp_bed,add_chr_prefix)
                self.make_picard_intervals(temp_bed,args.output,args.reference,add_chr_prefix)
                os.remove(temp_bed)
            else:
                self.make_picard_intervals(args.input,args.output,args.reference,add_chr_prefix) 