import os
import yaml
from pathlib import Path

# Get real paths for resources
RESOURCES_DIR=os.path.realpath(os.path.expanduser('~/resources'))
GENOMES_DIR=os.path.join(RESOURCES_DIR,'Genomes/Human/GRCh38')
INTERVALS_DIR=os.path.join(RESOURCES_DIR,'Interval_files')

# Default reference file
DEFAULT_REFERENCE=os.path.join(GENOMES_DIR,'Homo_sapiens.GRCh38.dna.primary_assembly.fa')

def get_reference_files(fasta_path):
    """Get reference files based on fasta path"""
    # Remove any compression extensions
    base=str(Path(fasta_path).with_suffix(''))
    # Handle both regular and gzipped files
    if base.endswith('.fa'):
        base=base[:-3]
    elif base.endswith('.fasta'):
        base=base[:-6]
    
    # Get dictionary file
    dict_file=f"{base}.dict"
    if not os.path.exists(dict_file):
        raise FileNotFoundError(f"Reference dictionary file not found: {dict_file}")
    
    # Get fai file
    fai_file=f"{fasta_path}.fai"
    if not os.path.exists(fai_file):
        raise FileNotFoundError(f"Reference index file not found: {fai_file}")
    
    return {
        'fasta':fasta_path,
        'dict':dict_file,
        'fai':fai_file
    }

# Standard reference paths
REFERENCE_PATHS={
    'GRCh38':{
        'fasta':DEFAULT_REFERENCE,
        'dict':os.path.join(GENOMES_DIR,'Homo_sapiens.GRCh38.dna.primary_assembly.dict'),
        'fai':os.path.join(GENOMES_DIR,'Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'),
        'chr_prefix':True
    },
    'hg38':{
        'fasta':DEFAULT_REFERENCE,
        'dict':os.path.join(GENOMES_DIR,'Homo_sapiens.GRCh38.dna.primary_assembly.dict'),
        'fai':os.path.join(GENOMES_DIR,'Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'),
        'chr_prefix':False
    }
}

class ConfigManager:
    def get_config_template(self,project_type,reference='GRCh38',reference_fasta=None):
        # Get reference files
        if reference_fasta:
            ref_files=get_reference_files(reference_fasta)
        else:
            ref_files={
                'fasta':REFERENCE_PATHS[reference]['fasta'],
                'dict':REFERENCE_PATHS[reference]['dict'],
                'fai':REFERENCE_PATHS[reference]['fai']
            }
        
        base_config={
            'project':{
                'name':'',
                'fastq_config':'fastq.yml',
                'bam_table':'bam.table',
                'sample_list':'sample.list'
            },
            'reference':{
                'fasta':ref_files['fasta'],
                'dict':ref_files['dict'],
                'fai':ref_files['fai'],
                'key':reference,
                'chr_prefix':REFERENCE_PATHS[reference]['chr_prefix']
            },
            'resources':{
                'library_key':'',
                'targets_key':'',
                'targets_bed':'',
                'targets_intervals':'',
                'picard_intervals':''
            }
        }
        
        if project_type=='somatic':
            base_config['project'].update({
                'tumor_list':'tumor.list',
                'normal_list':'normal.list',
                'pair_table':'pair.table'
            })
        
        return base_config
    
    def update_config(self,args):
        config=self.get_config_template(args.type,args.reference,args.reference_fasta)
        
        if args.project:
            config['project']['name']=args.project
        
        if args.targets:
            config['resources']['targets_key']=args.targets
            ref_key=args.reference.lower()
            config['resources']['targets_bed']=os.path.join(INTERVALS_DIR,f'{args.targets}-{ref_key}.bed')
            config['resources']['targets_intervals']=os.path.join(INTERVALS_DIR,f'{args.targets}-{ref_key}.intervals')
            config['resources']['picard_intervals']=os.path.join(INTERVALS_DIR,f'{args.targets}-{ref_key}.picard.intervals')
        
        with open(args.output,'w') as f:
            yaml.dump(config,f,default_flow_style=False) 