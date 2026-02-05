import pysam
import argparse
import csv
import sys
from collections import defaultdict

def parse_cigar_get_query_pos(cigar,ref_pos,read_start):
    """
    Parse CIGAR string to find query sequence position for a given reference position.
    
    Args:
        cigar: pysam CIGAR tuple
        ref_pos: Reference position (1-based)
        read_start: Read start position on reference (0-based)
    
    Returns:
        tuple: (query_pos: int or None, base_qual: int or None)
        query_pos is 0-based index into query sequence
    """
    if not cigar:
        return (None,None)
    
    ref_offset=ref_pos-1-read_start  # Convert to 0-based offset from read start
    query_pos=0
    
    for op,length in cigar:
        if op==0:  # M (match)
            if ref_offset<length:
                return (query_pos+ref_offset,None)
            query_pos+=length
            ref_offset-=length
        elif op==1:  # I (insertion)
            query_pos+=length
        elif op==2:  # D (deletion)
            if ref_offset<length:
                return (None,None)  # Position is in deletion
            ref_offset-=length
        elif op==3:  # N (skip)
            if ref_offset<length:
                return (None,None)  # Position is in skipped region
            ref_offset-=length
        elif op==4:  # S (soft clip)
            query_pos+=length
        elif op==5:  # H (hard clip)
            pass  # Doesn't consume query
        elif op==6:  # P (padding)
            pass
        elif op==7:  # = (match)
            if ref_offset<length:
                return (query_pos+ref_offset,None)
            query_pos+=length
            ref_offset-=length
        elif op==8:  # X (mismatch)
            if ref_offset<length:
                return (query_pos+ref_offset,None)
            query_pos+=length
            ref_offset-=length
        
        if ref_offset<0:
            break
    
    return (None,None)

def check_read_for_snvs(read,pos1,ref1,alt1,pos2,ref2,alt2):
    """
    Check if a read contains both SNVs at the expected positions.
    
    Args:
        read: pysam AlignedSegment
        pos1: First SNV position (1-based)
        ref1: Reference allele for SNV1
        alt1: Alternate allele for SNV1
        pos2: Second SNV position (1-based)
        ref2: Reference allele for SNV2
        alt2: Alternate allele for SNV2
    
    Returns:
        tuple: (has_snv1: bool, has_snv2: bool, snv1_base: str, snv2_base: str)
    """
    if not read.query_sequence:
        return (False,False,None,None)
    
    # Get query positions for both SNVs
    query_pos1,_=parse_cigar_get_query_pos(read.cigartuples,pos1,read.reference_start)
    query_pos2,_=parse_cigar_get_query_pos(read.cigartuples,pos2,read.reference_start)
    
    has_snv1=False
    has_snv2=False
    base1=None
    base2=None
    
    if query_pos1 is not None and query_pos1<len(read.query_sequence):
        base1=read.query_sequence[query_pos1].upper()
        if base1==alt1.upper():
            has_snv1=True
    
    if query_pos2 is not None and query_pos2<len(read.query_sequence):
        base2=read.query_sequence[query_pos2].upper()
        if base2==alt2.upper():
            has_snv2=True
    
    return (has_snv1,has_snv2,base1,base2)

def check_mnp_in_bam(bam_path,chr,pos1,ref1,alt1,pos2,ref2,alt2,min_mapq=0):
    """
    Check if two SNVs appear on the same reads in a BAM file.
    
    Args:
        bam_path: Path to BAM file
        chr: Chromosome name
        pos1: First SNV position (1-based)
        ref1: Reference allele for SNV1
        alt1: Alternate allele for SNV1
        pos2: Second SNV position (1-based)
        ref2: Reference allele for SNV2
        alt2: Alternate allele for SNV2
        min_mapq: Minimum mapping quality (default 0)
    
    Returns:
        dict: Statistics about read support
    """
    bam=pysam.AlignmentFile(bam_path,'rb')
    
    # Determine region to fetch (overlap both positions)
    start_pos=min(pos1,pos2)-1  # pysam uses 0-based, half-open
    end_pos=max(pos1,pos2)+1
    
    reads_with_snv1=0
    reads_with_snv2=0
    reads_with_both=0
    reads_with_neither=0
    total_reads=0
    
    read_details=[]
    
    try:
        for read in bam.fetch(chr,start_pos,end_pos):
            if read.mapping_quality<min_mapq:
                continue
            
            total_reads+=1
            has_snv1,has_snv2,base1,base2=check_read_for_snvs(read,pos1,ref1,alt1,pos2,ref2,alt2)
            
            if has_snv1:
                reads_with_snv1+=1
            if has_snv2:
                reads_with_snv2+=1
            if has_snv1 and has_snv2:
                reads_with_both+=1
            if not has_snv1 and not has_snv2:
                reads_with_neither+=1
            
            read_details.append({
                'read_name':read.query_name,
                'has_snv1':has_snv1,
                'has_snv2':has_snv2,
                'base1':base1,
                'base2':base2,
                'mapq':read.mapping_quality
            })
    finally:
        bam.close()
    
    return {
        'total_reads':total_reads,
        'reads_with_snv1':reads_with_snv1,
        'reads_with_snv2':reads_with_snv2,
        'reads_with_both':reads_with_both,
        'reads_with_neither':reads_with_neither,
        'cis_evidence':reads_with_both,
        'trans_evidence':reads_with_snv1+reads_with_snv2-2*reads_with_both,
        'read_details':read_details
    }

def process_mnp_pairs(mnp_csv,bam_path,output_csv,min_mapq=0):
    """
    Process MNP pairs from CSV and check read support in BAM.
    
    Args:
        mnp_csv: CSV file with MNP pairs from find_mnp_variants.py
        bam_path: Path to BAM file
        output_csv: Output CSV with read support information
        min_mapq: Minimum mapping quality
    """
    mnp_pairs=[]
    
    # Read MNP pairs from CSV
    with open(mnp_csv,'r') as infile:
        reader=csv.DictReader(infile)
        current_pair=None
        var1=None
        var2=None
        
        for row in reader:
            # Group pairs (they should be consecutive rows)
            if var1 is None:
                var1=row
            else:
                var2=row
                # Check if this is a pair (same sample, same gene, etc.)
                if (var1.get('Chr')==var2.get('Chr') and
                    var1.get('Sample.ID')==var2.get('Sample.ID') or
                    var1.get('Tumor.ID')==var2.get('Tumor.ID')):
                    mnp_pairs.append((var1,var2))
                var1=var2
                var2=None
    
    # Process each pair
    results=[]
    for i,(var1,var2) in enumerate(mnp_pairs):
        chr=var1.get('Chr','')
        pos1=int(var1.get('Start',0))
        ref1=var1.get('REF','')
        alt1=var1.get('ALT','')
        pos2=int(var2.get('Start',0))
        ref2=var2.get('REF','')
        alt2=var2.get('ALT','')
        
        print(f"Processing pair {i+1}/{len(mnp_pairs)}: {chr}:{pos1} {ref1}>{alt1} and {chr}:{pos2} {ref2}>{alt2}",file=sys.stderr)
        
        stats=check_mnp_in_bam(bam_path,chr,pos1,ref1,alt1,pos2,ref2,alt2,min_mapq)
        
        result={
            'Chr':chr,
            'Pos1':pos1,
            'Ref1':ref1,
            'Alt1':alt1,
            'Pos2':pos2,
            'Ref2':ref2,
            'Alt2':alt2,
            'Total_Reads':stats['total_reads'],
            'Reads_with_SNV1':stats['reads_with_snv1'],
            'Reads_with_SNV2':stats['reads_with_snv2'],
            'Reads_with_Both':stats['reads_with_both'],
            'Reads_with_Neither':stats['reads_with_neither'],
            'Cis_Evidence':stats['cis_evidence'],
            'Trans_Evidence':stats['trans_evidence'],
            'Cis_Ratio':f"{stats['reads_with_both']/max(stats['total_reads'],1):.3f}" if stats['total_reads']>0 else '0.000'
        }
        
        # Add original variant info
        for key in ['Gene','Sample.ID','Tumor.ID','HGVSc','HGVSp']:
            if key in var1:
                result[f'Var1_{key}']=var1.get(key,'.')
            if key in var2:
                result[f'Var2_{key}']=var2.get(key,'.')
        
        results.append(result)
    
    # Write output
    if output_csv:
        outfile=open(output_csv,'w')
    else:
        outfile=sys.stdout
    
    if results:
        fieldnames=results[0].keys()
        writer=csv.DictWriter(outfile,fieldnames=fieldnames,delimiter=',',
                             restval='.',extrasaction='ignore',quoting=csv.QUOTE_NONNUMERIC,
                             dialect='excel')
        writer.writeheader()
        for result in results:
            writer.writerow(result)
    
    if output_csv:
        outfile.close()
        print(f"Processed {len(results)} MNP pairs",file=sys.stderr)
        print(f"Output written to {output_csv}",file=sys.stderr)
    else:
        print(f"Processed {len(results)} MNP pairs",file=sys.stderr)

def main():
    parser=argparse.ArgumentParser(description='Check if MNP pairs appear on same reads in BAM file')
    parser.add_argument('-i','--mnp_csv',required=True,help='CSV file with MNP pairs from find_mnp_variants.py')
    parser.add_argument('-b','--bam',required=True,help='BAM file to check')
    parser.add_argument('-o','--output_csv',required=False,default=None,help='Output CSV with read support (default: stdout)')
    parser.add_argument('-q','--min_mapq',type=int,default=0,help='Minimum mapping quality (default: 0)')
    args=parser.parse_args()
    
    process_mnp_pairs(args.mnp_csv,args.bam,args.output_csv,args.min_mapq)

if __name__=='__main__':
    main()

