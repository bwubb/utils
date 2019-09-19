
import math
import csv
import sys

try:
    snakemake
except NameError:
    input=sys.argv[1:]
    output=[i.replace('.tumors.segments.filtered_CNV.features.tsv','.cns') for i in input if i.endswith('.tsv')]
else:
    input=snakemake.input
    output=snakemake.output


header=['chromosome','start','end','gene','log2','depth','weight','probes']
#lratio is weight
print(input)
print(output)
for i,file in enumerate(input):
    with open(input[i],'r') as infile, open(output[i],'w') as outfile:
        reader=csv.DictReader(infile,delimiter='\t',fieldnames=['chromosome','start','end','name','score','strand','gene'])
        writer=csv.DictWriter(outfile,delimiter='\t',fieldnames=header)
        writer.writeheader()
        outrow={}
        for row in reader:
            outrow['chromosome']=row['chromosome']
            outrow['start']=row['start']
            outrow['end']=row['end']
            outrow['gene']=row['gene']
            #info=row['name'].split(';')
            outrow.update(dict(i.split('=') for i in row['name'].split(';')))
            writer.writerow(outrow)