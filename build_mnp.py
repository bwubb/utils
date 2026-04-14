#!/usr/bin/env python3
import argparse
import subprocess
import sys
import vcfpy

def parse_id(vid):
    """Parse ID like chr22_16591614_C_T into (chr,pos,ref,alt)."""
    parts=vid.strip().split('_')
    if len(parts)<4:
        return None
    try:
        pos=int(parts[1])
        return (parts[0],pos,parts[2],parts[3])
    except (ValueError,IndexError):
        return None

def chrom_for_vcf(chr_):
    """Return contig name for VCF/ref (strip leading 'chr' for GRCh38-style)."""
    if chr_.lower().startswith('chr'):
        return chr_[3:]
    return chr_

def fetch_ref_region(ref_fasta,chr_,start,end):
    """Return reference sequence for chr:start-end (1-based inclusive)."""
    contig=chrom_for_vcf(chr_)
    region=f"{contig}:{start}-{end}"
    result=subprocess.run(
        ['samtools','faidx',ref_fasta,region],
        capture_output=True,text=True,check=True
    )
    seq=result.stdout.strip().split('\n')[-1]
    return seq

def make_header():
    """Minimal no-sample VCF header."""
    from vcfpy.header import Header,HeaderLine
    h=Header(lines=[HeaderLine('fileformat','VCFv4.3')],samples=None)
    return h

def load_pairs_from_input(input_path):
    """Read tab-separated pairs table (header: pair_id id1 id2 sample ...); return unique (id1,id2) pairs."""
    seen=set()
    with open(input_path) as f:
        for i,line in enumerate(f):
            toks=line.strip().split('\t')
            if len(toks)<3:
                continue
            if i==0 and toks[0]=='pair_id':
                continue
            id1,id2=toks[1],toks[2]
            if id1 and id2:
                seen.add((id1,id2))
    return list(seen)

def main():
    ap=argparse.ArgumentParser(description='Build no-sample MNP VCF from tested pairs table.')
    ap.add_argument('-r','--ref',required=True,help='Reference FASTA (indexed with samtools faidx)')
    ap.add_argument('-i','--input',required=True,help='Tab-separated pairs table: pair_id id1 id2 sample ...')
    ap.add_argument('-o','--output',required=True,help='Output VCF path')
    ap.add_argument('-H','--header',default=None,help='VCF file to load header from (optional; else minimal header)')
    args=ap.parse_args()

    if args.header:
        reader=vcfpy.Reader.from_path(args.header)
        header=reader.header.copy()
        reader.close()
    else:
        header=make_header()
    writer=vcfpy.Writer.from_path(args.output,header)
    n_adjacent=0
    n_2apart=0

    pairs=load_pairs_from_input(args.input)
    for id1,id2 in pairs:
        p1=parse_id(id1)
        p2=parse_id(id2)
        if not p1 or not p2:
            continue
        chr1,pos1,ref1,alt1=p1
        chr2,pos2,ref2,alt2=p2
        if chr1!=chr2:
            continue
        if pos1>pos2:
            pos1,pos2=pos2,pos1
            ref1,ref2=ref2,ref1
            alt1,alt2=alt2,alt1

        if pos2==pos1+1:
            ref=ref1+ref2
            alt=alt1+alt2
            pos=pos1
            n_adjacent+=1
        elif pos2==pos1+2:
            seq=fetch_ref_region(args.ref,chr1,pos1,pos2)
            if len(seq)!=3:
                continue
            if seq[0]!=ref1 or seq[2]!=ref2:
                continue
            ref=seq
            alt=alt1+seq[1]+alt2
            pos=pos1
            n_2apart+=1
        else:
            continue

        vid=f"{chr1}_{pos1}_{ref}_{alt}"
        alt_rec=vcfpy.Substitution(type_='MNV',value=alt)
        rec=vcfpy.Record(
            CHROM=chrom_for_vcf(chr1),
            POS=pos,
            ID=[vid],
            REF=ref,
            ALT=[alt_rec],
            QUAL=None,
            FILTER=['PASS'],
            INFO={},
            FORMAT=None,
            calls=None
        )
        writer.write_record(rec)

    writer.close()
    print(f"Wrote {n_adjacent} adjacent + {n_2apart} 2-apart MNPs -> {args.output}",file=sys.stderr)

if __name__=='__main__':
    main()
