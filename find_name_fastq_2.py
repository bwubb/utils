#!/usr/bin/env python

import argparse
import glob
import csv
import os


def get_args():
	p = argparse.ArgumentParser()
	p.add_argument('-i', '--infile', help='Sequence Project submission file')
	p.add_argument('-f', '--fastq_dir', help='Directory where fastq files are located')
	p.add_argument('-d', '--dest_dir', help='Destination directory')
	p.add_argument('-R', '--rename', action='store_true', default=False, help='Rename files. Default(False)')
	args = p.parse_args()
	print("Arguments selected:")
	for a, b in vars(args).items():
		print('{}: {}'.format(a,b))
	return args

def find_fastqs(fqdir, run, lane, index1,index2):
	fqs = sorted(glob.glob('{0}/{1}*_{2}_[1-2]_{3}-{4}.fastq.gz'.format(fqdir,run,lane,index1,index2)))
	if len(fqs) != 2:
		print("Error:", fqs, "found matching pattern",run,lane,index1,index2)
		return None
	else:
		return fqs

def new_name(fq, sample):
	x = os.path.basename(fq).split('_')
	if x[1] == 's':
		x.pop(1)
	x.insert(0,sample)
	return '_'.join(x)

def work(file, fqdir, dest, rename):
	m = file.split('.')
	m.insert(-1, 'missing')
	with open(file, 'r') as infile, open('.'.join(m), 'w') as missing:
		reader = csv.reader(infile, delimiter='\t')
		writer = csv.writer(missing, delimiter='\t')
		for row in reader:
			if row[0].startswith('Run'):
				writer.writerow(row)
				continue
			fqs = find_fastqs(fqdir, *row[:4])
			if not fqs:
				writer.writerow(row)
				continue
			for fq in fqs:
				print(fq, '>', '/'.join([dest, new_name(fq, row[4])]))
				if rename:
					os.rename(fq, '/'.join([dest, new_name(fq, row[4])]))
			
def main(argv=None):
	args = get_args()
	work(os.path.abspath(args.infile), os.path.abspath(args.fastq_dir), os.path.abspath(args.dest_dir), args.rename)


if __name__ == '__main__':
	main()