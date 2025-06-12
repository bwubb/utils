#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from .samples import SampleManager
from .config import ConfigManager
from .intervals import IntervalManager

def get_args(argv=None):
    p=argparse.ArgumentParser()
    subparsers=p.add_subparsers(dest='command',help='Available commands')
    
    # Update samples
    samples=subparsers.add_parser('update_samples',help='Update sample lists and tables')
    samples.add_argument('-i','--input',help='Input FASTQ directory or sample submission file')
    samples.add_argument('-o','--output',help='Output directory for sample files')
    samples.add_argument('--source',choices=['fgc','nextseq'],default='fgc',help='Source of FASTQ files')
    samples.add_argument('--action',choices=['copy','rename','dryrun'],default='dryrun',help='Action to take')
    samples.add_argument('--type',choices=['germline','somatic'],required=True,help='Project type')
    
    # Update config
    config=subparsers.add_parser('update_config',help='Update project configuration')
    config.add_argument('-i','--input',help='Input config file')
    config.add_argument('-o','--output',help='Output config file')
    config.add_argument('--project',help='Project name')
    config.add_argument('--reference',choices=['GRCh38','hg38'],default='GRCh38',help='Reference genome')
    config.add_argument('--targets',help='Target regions')
    config.add_argument('--type',choices=['germline','somatic'],required=True,help='Project type')
    
    # Process intervals
    intervals=subparsers.add_parser('process_intervals',help='Process interval files')
    intervals.add_argument('-i','--input',required=True,help='Input file (BED or intervals)')
    intervals.add_argument('-o','--output',required=True,help='Output file')
    intervals.add_argument('-r','--reference',choices=['GRCh38','hg38'],default='GRCh38',help='Reference genome')
    intervals.add_argument('--format',choices=['bed','intervals','picard'],required=True,help='Output format')
    
    return p.parse_args(argv)

def main(argv=None):
    args=get_args(argv)
    
    if args.command=='update_samples':
        manager=SampleManager()
        manager.update_samples(args)
    
    elif args.command=='update_config':
        manager=ConfigManager()
        manager.update_config(args)
    
    elif args.command=='process_intervals':
        manager=IntervalManager()
        manager.process_intervals(args)

if __name__=='__main__':
    main() 