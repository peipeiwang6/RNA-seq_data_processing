import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for RNA-seq data processing')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-SRA', help='SRA sample ID', required=True)
	req_group.add_argument('-genome_seq', help='genome sequence fasta file', required=True)
	req_group.add_argument('-gff', help='gff3 file', required=True)
	req_group.add_argument('-layout', help='PE: paired end, SE: single end', required=True)
	req_group.add_argument('-workdir', help='the path to your workdir', required=True)
	
	# optional
	inp_group = parser.add_argument_group(title='OPTIONAL INPUT')
	inp_group.add_argument('-trim', help='trim the reads or not, y or n',default='n')
	inp_group.add_argument('-threads', help='how many threads to be used',default=4)
	inp_group.add_argument('-adapters', help='adapter sequences to be trimmed',default='None')
	inp_group.add_argument('-seedMis', help='seed mismatches: specifies the maximum mismatch count which will still allow a full match to be performed',default=2)
	inp_group.add_argument('-pClipThres', help='palindromeClipThreshold: specifies how accurate the match between the two adapter ligated reads must be for PE palindrome read alignment.',default=30)
	inp_group.add_argument('-sClipThres', help='simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.',default=10)
	inp_group.add_argument('-LEADING', help='Specifies the minimum quality required to keep a base',default=3)
	inp_group.add_argument('-TRAILING', help='Specifies the minimum quality required to keep a base.',default=3)
	inp_group.add_argument('-windowSize', help='specifies the number of bases to average across',default=4)
	inp_group.add_argument('-requiredQuality', help='specifies the average quality required',default=20)
	inp_group.add_argument('-MINLEN', help='Specifies the minimum length of reads to be kept.',default=36)
	inp_group.add_argument('-time', help='hours asked for job running',default=4)
	inp_group.add_argument('-ntasks', help='ntasks asked for job running',default=1)
	inp_group.add_argument('-cpus', help='cpus-per-task asked for job running',default=8)
	inp_group.add_argument('-mem', help='mem asked for job running',default=20)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	
	
	SRA = args.SRA
	slum_code = open('Hisat_%s.sh'%SRA, 'w')
	slum_code.write('#!/bin/sh --login\n#SBATCH --time=%s:00:00\n#SBATCH --ntasks=%s\n#SBATCH --cpus-per-task=%s\n#SBATCH --mem=%sG\n#SBATCH --job-name Hisat_%s.sh\n#SBATCH -e Hisat_%s.sh.e\n#SBATCH -o Hisat_%s.sh.o\ncd %s\n'%(args.time,args.ntasks,args.cpus,args.mem,SRA,SRA,SRA,args.workdir))
	slum_code.write('module load SRA-Toolkit FastQC Trimmomatic hisat2 SAMtools GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA HTSeq\n')
	
#	# build the index for the genome
#	slum_code.write('hisat2-build %s %s_gi\n'%(args.genome_seq,args.genome_seq))
		
	if args.layout == 'PE':
		slum_code.write('fastq-dump --split-files %s\n'%SRA)
		if args.trim == 'y':
			# check the quality before trimming
			slum_code.write('fastqc -f fastq %s_1.fastq\n'%SRA)
			slum_code.write('fastqc -f fastq %s_2.fastq\n'%SRA)
			slum_code.write('java -jar /opt/software/Trimmomatic/0.36-Java-1.8.0_92/trimmomatic-0.36.jar PE -threads %s %s_1.fastq %s_2.fastq %s_1.trimP %s_1.trimU %s_2.trimP %s_2.trimU ILLUMINACLIP:%s:%s:%s:%s LEADING:%s TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s\n'%(args.threads, SRA, SRA, SRA, SRA, SRA, SRA, args.adapters, args.seedMis, args.pClipThres, args.sClipThres, args.LEADING,args.TRAILING, args.windowSize, args.requiredQuality, args.MINLEN))
			# check the quality after trimming
			slum_code.write('fastqc -f fastq %s_1.trimP\n'%SRA)
			slum_code.write('fastqc -f fastq %s_2.trimP\n'%SRA)
			# map the reads to the genome
			slum_code.write('hisat2 -p 4 --dta -x %s_gi -1 %s_1.trimP -2 %s_2.trimP -S %s.sam\n'%(args.genome_seq, SRA, SRA, SRA))
		else:
			# map the reads to the genome
			slum_code.write('hisat2 -p 4 --dta -x %s_gi -1 %s_1.trimP -2 %s_2.trimP -S %s.sam\n'%(args.genome_seq, SRA, SRA, SRA))

	if args.layout == 'SE':
		slum_code.write('fastq-dump %s\n'%SRA)
		if args.trim == 'y':
			# check the quality before trimming
			slum_code.write('fastqc -f fastq %s.fastq\n'%SRA)
			slum_code.write('java -jar /opt/software/Trimmomatic/0.36-Java-1.8.0_92/trimmomatic-0.36.jar SE -threads %s %s.fastq %s.trimP %s.trimU ILLUMINACLIP:%s:%s:%s:%s LEADING:%s TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s\n'%(args.threads, SRA, SRA, SRA, args.adapters, args.seedMis, args.pClipThres, args.sClipThres, args.LEADING,args.TRAILING, args.windowSize, args.requiredQuality, args.MINLEN))
			# check the quality after trimming
			slum_code.write('fastqc -f fastq %s.trimP\n'%SRA)
			# map the reads to the genome
			slum_code.write('hisat2 -p 4 --dta -x %s_gi -U %s.trimP -S %s.sam\n'%(args.genome_seq, SRA, SRA))
		else:
			# map the reads to the genome
			slum_code.write('hisat2 -p 4 --dta -x %s_gi -U %s.trimP -S %s.sam\n'%(args.genome_seq, SRA, SRA))
		
	# sort the sam file
	slum_code.write('samtools sort  -n -O sam, %s.sam  -o %s_sorted.sam\n'%(SRA, SRA))
	# get uniquely mapped reads 
	slum_code.write('python 02_keep_reads_with_quality_60_and_unique_mapping.py %s_sorted.sam\n'%(SRA))
	# get read counts
	slum_code.write('htseq-count --format=sam -m union -s no -t gene -i ID %s_sorted_quality_60_unique.sam %s > HTSeqCount_%s.out\n'%(SRA, args.gff, SRA))
	
	slum_code.close()

if __name__ == '__main__':
	main()
	
	
