import os
import sys
import csv
from collections import defaultdict

#CODEX2

#Need input check
try:
    snakemake
except NameError:
    input=sys.argv[1:]
    output=[i.replace('.txt','_CNV.bed') for i in input if i.endswith('.txt')]
else:
    input=snakemake.input
    output=snakemake.output

bed_fields=['chrom','chromStart','chromEnd','name','score','strand']
chr_sort={**{str(c):c for c in list(range(1,23))},**{'X':23,'Y':24}}
segments=defaultdict(list)
for i,file in enumerate(input):
    with open(file,'r') as segments_file:
        print(f'Reading {segments_file.name}...')
        reader=csv.DictReader(segments_file,delimiter='\t')
        for row in reader:
            segments[row['sample_name']].append({key:row[key] for key in ['chr','st_bp','ed_bp','cnv','copy_no']})
    samples=segments.keys()
    #
    for sample in samples:
        sample_segs=sorted(segments[sample],key=lambda x: (chr_sort[x['chr']],int(x['st_bp'])))
        try:
            snakemake
        except NameError:
            out_name=f'{os.path.dirname(output[i])}/{sample}.codex2.{os.path.basename(output[i])}'
        else:
            out_name=(output[f] for f in output if os.path.basename(f)==f'{sample}.codex2.segments_CNV.bed')
        with open(out_name,'w') as bed_file:
            #print('Writing')
            writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields,lineterminator='\n')
            for segment in sample_segs:
                start=int(segment['st_bp'])-1
                end=int(segment['ed_bp'])
                if start>=end:
                    start=int(segment['ed_bp'])-1
                    end=int(segment['st_bp'])
                    strand='-'
                else:
                    strand="+"
                try:
                    length=end-start+1
                except ValueError:
                    length='.'
                name=f"{length}bp;{segment['cnv']}"#consider adding gain,amp,del,loss
                score=f"{float(segment['copy_no'])*100:.{0}f}"
                bed_row={'chrom':segment['chr'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':f"{name}",'score':f"{score}",'strand':strand}
                writer.writerow(bed_row)


