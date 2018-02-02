
import glob
import os

fastqs=glob.glob('FASTQ/*.fastq.gz')
samples=set()
for file in fastqs:
    if not os.path.basename(file).startswith('FGC'):
        samples.add(os.path.basename(file).split('_')[0])

with open('samples.list','w') as file:
    file.writelines('%s\n' % x for x in sorted(list(samples)))
        