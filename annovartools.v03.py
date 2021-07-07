import argparse
import vcfpy
import csv
import os
import datetime
import yaml
import re
import itertools
import sys
from collections import defaultdict


'''
Brad Wubbenhorst
bwubb@pennmedicine.upenn.edu
https://github.com/bwubb
'''

class TranscriptError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class Func_refGeneError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class SampleError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class Multianno(object):
    def __init__(self,record,transcript,header=None):
        '''
        record = vcfpy.Record
        record.INFO = vcfpy.OrderedDict that contains annotations
        transcript = {k:[v]}
        '''
        #######################################################################
        #TO DO:
        #ncRNA_regions
        #Improve intronic annotation? HGVS numbering from refGene features | Should be ANNOVAR feature...
        #Anti-target region annotation
        #######################################################################
        self.fields=defaultdict(list)
        self._pfields=defaultdict(list)
        self.fields['Chr']=record.CHROM
        self.fields['Start']=str(record.POS)
        self.fields['Ref']=record.REF
        self.fields['Alt']=[x.value for x in record.ALT]#vcfpy.AltRecord
        self.fields['ID']=record.ID
        self.fields['FILTER']='|'.join(record.FILTER)
        self.fields['Gene.refGene']=record.INFO['Gene.refGene']
        Func_refGene={'exonic':self.exonic_func,
                'splicing':self.splicing_func,
                'exonic\\3xbsplicing':self.exonicsplicing_func,
                'UTR3':self.utr_func,
                'UTR5':self.utr_func,
                'introic':self.intronic_func,
                'introic\\3xbsplicing':self.intronic_func
                #'intergenic':intergenic_func,
                }
        #Args_refGene dict to pass specific required annotation values?
        for key in record.INFO['Func.refGene']:
            self.fields['GenomicRegion.refGene'].append(key)
            Func_refGene.get(key,self.UNKNOWN_func)(record.INFO,transcript)
        
        ### Step 2
        self.other_fields(record.INFO)
        #TO DO: Calculate END for SV feature
        
        ### Step 3
        #For VCF only
        self.call_record={}
        for call in record.calls:
            self.call_record[call.sample]=self.call_fields(call)
        self.sample=record.calls[0].sample
        #What about no sample vcf?
    
    def UNKNOWN_func(self,INFO,transcript):
        pass
    
    @staticmethod
    def drop_gene(Annot,GENE_KEY):
        rep=dict((re.escape(f'{gene}:'),'') for gene in GENE_KEY)
        pattern=re.compile("|".join(f'^{k}' for k in rep.keys()))
        return pattern.sub(lambda m: rep[re.escape(m.group(0))],Annot)
    
    def exonic_func(self,INFO,transcript):
        def exonic_update(fields,vals):
            try:#Nonframeshift sub have no aa, <INSERT> custom fix
                fields['AAChange.refGene'].append(vals[3])
            except IndexError:
                fields['AAChange.refGene'].append('.')
            try:
                fields['NTChange.refGene'].append(vals[2])
                fields['Exon.refGene'].append(vals[1].replace('exon','Exon '))
            except IndexError:
                if vals[1]=='wholegene':
                    fields['NTChange.refGene'].append('.')
                    fields['Exon.refGene'].append(vals[1])
            if vals[0].startswith('NM'):
                fields['Transcript.refGene'].append(vals[0])
            return fields
        annotation=INFO.get('AAChange.refGene',list())
        
        for Annot in annotation:
            gene_detail=self.drop_gene(Annot,transcript.keys())
            self.fields['GeneDetail.refGene'].append(gene_detail)
            vals=gene_detail.split(':')
            self.fields=exonic_update(self.fields,vals)
            if self.check_set(vals,itertools.chain.from_iterable(transcript.values())):
                self._pfields=exonic_update(self._pfields,vals)
    
    def splicing_func(self,INFO,transcript):
        def splicing_update(fields,vals):
            fields['NTChange.refGene'].append(vals[2])
            fields['Exon.refGene'].append(vals[1].replace('exon','Exon '))
            fields['Transcript.refGene'].append(vals[0])
            return fields
        
        annotation=self.split_x3b(INFO.get('GeneDetail.refGene',list()))
        for Annot in annotation:
            gene_detail=self.drop_gene(Annot,transcript.keys())
            self.fields['GeneDetail.refGene'].append(gene_detail)
            vals=gene_detail.split(':')
            self.fields=splicing_update(self.fields,vals)
            if self.check_set(vals,itertools.chain.from_iterable(transcript.values())):
                self._pfields=splicing_update(self._pfields,vals)
    
    def exonicsplicing_func(self,INFO,transcript):
        self.exonic_func(INFO,transcript)
        self.splicing_func(INFO,transcript)
    
    def intronic_func(self,INFO,transcript):
        #intronic is no longer annovar'd well. No transcript given.
        #INFO could equal dict of GeneDetail + dbscSNV
        for x in ['dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE']:
            try:
                float(INFO.get(x))
            except ValueError:
                assert INFO.get(x,['.']) in [[''],['.'],list()]
                pass
            else:
                #print('intronic\\x3bsplicing')
                self.fields['GenomicRegion.refGene'][-1]='intronic\\x3bsplicing'
    
    def utr_func(self,INFO,transcript):
        def UTR_update(fields,vals):
            fields['NTChange.refGene'].append(vals[1])
            fields['Transcript.refGene'].append(vals[0])
            return fields
        
        annotation=self.split_x3b(INFO.get('GeneDetail.refGene',list()))
        for Annot in annotation:
            gene_detail=self.drop_gene(Annot,transcript.keys())
            self.fields['GeneDetail.refGene'].append(gene_detail)
            vals=gene_detail.split(':')
            self.fields=UTR_update(self.fields,vals)
            if self.check_set(vals,itertools.chain.from_iterable(transcript.values())):
                self._pfields=UTR_update(self._pfields,vals)
    
    #ncRNA_exonic Need for TMB.
    
    def intergenic_func(self,INFO,transcript):
        #probably upstream and downstream too
        #   #   #GeneDetail.refGene
        #intergenic	APC;SRP19	dist=2211;dist=12738
        #intergenic	APC;SRP19	dist=6158;dist=8791
        pass
    
    def other_fields(self,INFO,header=None):
        '''
        For headers/keys within INFO that are not refGene annotations
        '''
        default=['SampleID','TumorID','NormalID','CHR','START','REF','ALT','FILTER']
        refGene=[k for k in self.fields.keys() if k.endswith('refGene')]
        #parse header to catch db versions
        #if not in header; pass
        for key,val in INFO.items():
            if any([(key in default),(key in refGene)]):
                continue
            else:
                #commans in FILTER slipped through, This is a temp fix for a larger issue
                #Where are the values coming from that I want to display.
                self.fields[key]=val
    
    @staticmethod
    def call_fields(call):
        '''
        VCF FORMAT Fields can no longer be trusted;
        Certain softwares think they are above convention
        This will attempt to format [Zyg, GT, DP, AAD, AAF]
        '''
        def freq_thres(f):
            #not called trigger from warnings/exceptions
            #additional dict k/v
            return float(f)>0.03
        
        def Zyg(GT):
            '''
            Assign Zygosity designation based on observed alleles
            '''
            if GT==['.','.']:
                return '.'
            elif sum(GT)==0:
                return 'HOM_REF'
            elif len(set(GT))==1:
                return 'HOM_ALT'
            elif len(set(GT))==2:
                return 'HET_ALT'
            elif len(set(GT))>2:
                return 'COMPLEX'
            else:
                #need actual undetermined method
                return '.'
        
        def CSV(l):
            return ';'.join(map(str,l))
        
        try:
            DP=int(call.data.get('DP',0))
        except TypeError:
            DP=0
        try:
            AD=call.data.get('AD',[0,0])
        except TypeError:
            AD=[0,0]
        #varscan2 fix
        #if int(call.data.get('RD',-1))>=0:
        #    AD=[call.data['RD'],AD]
        #warnings exceptions over missing values
        #This has gotten so convoluted I need to find a better way of going through this.
        #I need a bcftools method for fixing FORMAT and VALUES
        if not all([type(AD)==list,len(AD)>1]):
            AD=[0,0]
        if DP==0:
            try:
                DP=sum(AD)
            except TypeError as e:
                print(f"{e}: AD = {call.data.get('AD',[0,0])}")
        try:
            AF=[f'{x/DP:.3f}' for x in map(float,AD)]
        except (ZeroDivisionError,TypeError) as e:
            #print(f'{e}: DP = {str(DP)}, AD = {str(AD)}')
            DP=0
            AF=['0.0' for x in range(0,len(AD))]
        
        GT=[]
        for i,freq in enumerate(AF):
            if freq_thres(freq):
                GT.append(i)
        if GT==[]:
            GT.append('.')
        if len(GT)==1:
            GT.append(GT[0])
        
        #Need functionality for ./.
        #   had {True:i,False:'.'} but then need to assign GTs
        #def return ','.join(map(str,x))
        return [Zyg(GT),CSV(GT),str(DP),CSV(AD[1:]),CSV(AF[1:])]
    
    @staticmethod
    def check_set(test,acceptable):
        if len(set(test).intersection(set(acceptable)))>0:
            return True
        return False
    
    @staticmethod
    def split_x3b(details):
        proper=[]
        for item in details:
            assert type(item)==str
            for x in item.split('\\x3b'):
                proper.append(x)
        return proper
    
    def unhex_format_fields(self):
        '''
        Turns lists into strings
        Replaces ('\\x3d','=') ('\\x3b',';')
        Converts empty string to .
        returns dict()
        '''
        def string_replace(x):
            v=str(x)
            if type(x)==list:
                v='|'.join(map(str,x))#
            v=v.replace('\\x3d','=').replace('\\x3b',';')
            if v=='':
                v='.'
            return v
        nohex={}
        for key,val in self.fields.items():
            nohex[key]=string_replace(val)
        if hasattr(self,'_pfields'):
            for key,val in self._pfields.items():
                nohex[key]=string_replace(val)
        return nohex #This is a dictionary
    
    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return str(self.unhex_format_fields())
    
    def __call__(self,sample=None,sample2=None):
        '''
        returns dict unhexed,formatted values
        for single sample or tumor,normal
        '''
        if sample==None:
            sample=self.sample
        fields=self.unhex_format_fields()
        if sample2:
            #paired
            h1=['Tumor_Zyg','Tumor_Genotype','Tumor_Total_Depth','Tumor_ALT_AlleleDepth','Tumor_ALT_AlleleFrac']
            if 'TUMOR' in self.call_record.keys():
                t=self.call_record['TUMOR']
            else:
                t=self.call_record[sample]
            fields.update(dict(zip(h1,t)))
            fields['TumorID']=sample
            h2=['Normal_Zyg','Normal_Genotype','Normal_Total_Depth','Normal_ALT_AlleleDepth','Normal_ALT_AlleleFrac']
            if 'NORMAL' in self.call_record.keys():
                n=self.call_record['NORMAL']
            else:
                n=self.call_record[sample2]
            fields.update(dict(zip(h2,n)))
            fields['NormalID']=sample2
        elif sample:
            #single
            h=['Zyg','Genotype','Total_Depth','ALT_AlleleDepth','ALT_AlleleFrac']
            i=self.call_record[sample]
            fields.update(dict(zip(h,i)))
            fields['SampleID']=sample
        else:
            pass
        return fields

### Functions ###

def report_items(d):
    for k,v in d.items():
        print(k,':',v)

def load_TRANSCRIPT(TRANSCRIPT):
    with open(TRANSCRIPT,'r') as load_ptx:
        print(f'Loading {load_ptx.name}')
        #If yaml do this; else
        pTX=yaml.load(load_ptx,Loader=yaml.BaseLoader)
    return pTX

def parse_arguments():
    p=argparse.ArgumentParser()
    p.add_argument('-I','--input_fp',help='Annovar multianno.vcf')
    p.add_argument('-O','--output_fp',help='Output report.tsv')#Can default to sys.stdout which has .write()
    p.add_argument('-m','--mode',choices=['single','paired','cohort'],default='single')
    p.add_argument('--header',default='annovar_header.txt',help='annovar header file')
    p.add_argument('--transcript_yaml',default='/home/bwubb/resources/Transcript_files/refGene_transcripts.20180426.yaml',help='Principal transcipts in yaml form.')
    #p.add_argument('--all-zyg',default=False,action='store_true',help='Turn off HOM_REF variant filter. Default behavior removes HOM_REF variants.')
    p.add_argument('samples',nargs=argparse.REMAINDER,help="Optional subset of samples. If 'paired' mode, use 'Tumor' 'Normal'.")
    return vars(p.parse_args())

def parse_snakemake():
    argv={}
    argv['input_fp']=snakemake.input[0]
    argv['output_fp']=snakemake.output[0]
    argv['mode']=snakemake.params['mode']
    argv['transcript_yaml']=snamemake.params['transcript_yaml']
    argv['header']=snakemake.params['header']
    argv['all_']
    argv['samples']=snakemake.params['samples']
    print(argv)
    return argv

def default_data():
    #load default
    header=['Chr','Start','Ref','Alt']
    header+=['GenomicRegion.refGene','ExonicFunc.refGene','Gene.refGene','Transcript.refGene']
    header+=['Exon.refGene','NTChange.refGene','AAChange.refGene']
    header+=['Zyg','Genotype','Total_Depth','ALT_AlleleDepth','ALT_AlleleFrac']
    INPUT='test.somatic.pass.hg19_multianno.vcf'
    OUTPUT='test.somatic.pass.hg19_multianno.report.tsv'

### ### ### ###

def main(argv=None):
    print('Arguments Initialized...')
    for k,v in argv.items():
        print(k,':',v)
    INPUT=os.path.abspath(argv['input_fp'])
    OUTPUT=open(os.path.abspath(argv['output_fp']),'w')
    with open(argv['header'],'r') as file:
        header=file.read().splitlines()
    
    ANNOTATIONS=defaultdict(int)
    TRANSCRIPT=argv['transcript_yaml']
    pTX=load_TRANSCRIPT(TRANSCRIPT)
    
    writer=csv.DictWriter(OUTPUT,delimiter='\t',fieldnames=header,restval='.',extrasaction='ignore')
    writer.writeheader()
    print(header)
    VcfReader=vcfpy.Reader.from_path(INPUT)
    for record in VcfReader:
        transcript={}
        for gene in record.INFO['Gene.refGene']:
            transcript[gene]=[pTX.get(gene,'')]
        #print(transcript)
        if 'UNKNOWN' in record.INFO.get('AAChange.refGene'):
            #Hard Skip1
            continue
        DATA=Multianno(record,transcript)
        if pTX.get(gene,'')!='' and len(set(record.INFO['Gene.refGene']))==1 and pTX.get(gene,'') not in DATA.fields['Transcript.refGene']:
            #Hard skip2
            #It was found annoying to have transctipt1|transcript2 and niether of them were the princpal one.
            continue
        if argv['mode']=='single':
            if argv['samples']==None:
                writer.writerow(DATA())
                ANNOTATIONS[DATA.fields['GenomicRegion.refGene'][0]]+=1
            else:
                for sample in argv['samples']:
                    #try except SampleError
                    writer.writerow(DATA(sample))
                    ANNOTATIONS[DATA.fields['GenomicRegion.refGene'][0]]+=1
        elif argv['mode']=='paired':
            tumor=argv['samples'][0]
            normal=argv['samples'][1]
            #debug ^
            writer.writerow(DATA(tumor,normal))
            ANNOTATIONS[DATA.fields['GenomicRegion.refGene'][0]]+=1
        #Warning about fields not found in header
        elif argv['mode']=='cohort':
            called_samples=[call.sample for call in record.calls if call.called]
            c=0
            for sample in called_samples:
                if DATA.call_record[sample][0]!='HOM_REF':
                    #call_record is stored as list
                    c+=1
                    writer.writerow(DATA(sample))
            print(f'{c} samples')
            ANNOTATIONS[DATA.fields['GenomicRegion.refGene'][0]]+=1
        #elif argv['mode']=='no_sample':
            
    OUTPUT.close()
    report_items(ANNOTATIONS)
    


if __name__=='__main__':
    csv.field_size_limit(sys.maxsize)
    main(parse_arguments())
    #try:
    #    snakemake
    #except NameError:
    #    main(parse_arguments())
    #else:
    #    main(parse_snakemake())
