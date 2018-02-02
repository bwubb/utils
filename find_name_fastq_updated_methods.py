


import csv
import re
import os
from collections import defaultdict

workdir='/project/ngsc_data/PI_INVESTIGATIONS/Nathanson_Katherine'
project='KNKD-ChaoWES'

aaa-studyinfo=defaultdict(list)
with open(os.path.join(workdir,project,'AAA-StudyInfo.txt'),'r') as AAA:
    reader=csv.DictReader(AAA,delimiter='\t')
    for row in reader:
        aaa-studyinfo[row['SAMP_name']].append({'run':row['RULA_Run'],'lane':row['RULA_Lane'],'index':row['RULA_Barcode'],'status':row['RULA_Status']})
    


    #{
    #"sample1":[{
    #"run":
    #"lane":
    #"index":str().replace(',','-')
    #"status":
    #},{}]
    #}

def fastq_from_dir(dir='./FASTQ',verbose=False):
    fastq_info=defaultdict(list)
    prog=re.compile('([\w\-\+]+)_(FGC\d{4})_([1-8]{1})_1_([ACGT\-]+)\.fastq\.gz')
    for file in os.listdir(dir):
        fastq=prog.match(file)
        if fastq:
            fastq_info[fastq.group(1)].append({'run':fastq.group(2),'lane':fastq.group(3),'index':fastq.group(4),'status':'dummy_status'})
        else:
            if verbose:
                print(file,'did not match.')
    return fastq_info