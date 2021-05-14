import sys,os
for files in os.listdir('./'):
	if files.endswith('_1.fastq.trimP'):
		name = files.split('_1.fastq.trimP')[0]
		out = open('Hisat2_%s.sh'%name,'w')
		out.write('#!/bin/sh --login\n#SBATCH --time=4:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=20G\n#SBATCH --job-name Hisat2_%s.sh\n#SBATCH -e Hisat2_%s.sh.e\n#SBATCH -o Hisat2_%s.sh.o\ncd /mnt/scratch/peipeiw/Allelic_specific/05_hisat2/Other_SRAs\n'%(name,name,name))
		out.write('module load SRA-Toolkit\nmodule load FastQC\nmodule load hisat2\nmodule load SAMtools\n')
		# hisat2
		out.write("hisat2 -p 4 --dta -x /mnt/scratch/peipeiw/Allelic_specific/05_hisat2/Pvirgatum_516_v5.0_gi -1 %s_1.fastq.trimP -2 %s_2.fastq.trimP -S %s.sam\n" % (name,name,name))
		# sort the sam file
		out.write("samtools sort  -n -O sam, %s.sam  -o %s_sorted.sam\n"%(name,name))
		# htseq, read counts
		out.write('module load GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA\nmodule load HTSeq\n')
		out.write('htseq-count --format=sam -m union -s no -t gene -i ID %s_sorted.sam /mnt/scratch/peipeiw/Allelic_specific/05_hisat2/Pvirgatum_516_v5.1.gene.gff3 > HTSeqCount_%s.out\n'%(name,name))
		out.close()

