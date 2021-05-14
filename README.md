# RNA-seq_data_processing

## you can apply the code 01_RNA-seq_processing.py to do the jobs automatically, or go through each step one by one. An example command line for this code is as follow, which will produce a .sh file to be submitted to the slurm queue

> python 01_RNA_seq_processing.py -SRA SRR14066076 -genome_seq Pvirgatum_516_v5.0.fa -gff Pvirgatum_516_v5.1.gene.gff3 -layout PE -workdir /mnt/scratch/peipeiw/Data_for_Kenia/For_RNA_seq_pipeline -trim y -adapters all_PE_adapters.fa


## load modules
> module load SRA-Toolkit/2.10.7-centos_linux64
### save the config for sra toolkit, otherwise the fastq-dump won't work
> vdb-config --interactive  
> module load FastQC
> module load Trimmomatic
> module load hisat2
> module load SAMtools
> module load GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA
> module load HTSeq

## download and split the fastq file from SRA
> fastq-dump --split-files SRR14066076

## check the quality, QC
> fastqc -f fastq SRR14066076_1.fastq
> fastqc -f fastq SRR14066076_2.fastq

## trim the fastq file
> java -jar /opt/software/Trimmomatic/0.36-Java-1.8.0_92/trimmomatic-0.36.jar PE -threads 4 SRR14066076_1.fastq SRR14066076_2.fastq SRR14066076_1.trimP SRR14066076_1.trimU SRR14066076_2.trimP SRR14066076_2.trimU ILLUMINACLIP:all_PE_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

## build the genome index
> hisat2-build YourGenome.fa YourGenome_gi

## map the read
> hisat2 -p 4 --dta -x YourGenome_gi -1 SRR14066076_1.fastq -2 SRR14066076_2.fastq -S SRR14066076.sam

## sort the sam file
> samtools sort  -n -O sam, SRR14066076.sam  -o SRR14066076_sorted.sam

## get count
> htseq-count --format=sam -m union -s no -t gene -i ID SRR14066076_sorted.sam YourGenome.gff3 > HTSeqCount_SRR14066076.out

## combine the read count files
> python 08_combine_read_counts.py

## call the TPM, which will be done in R
> Rscript 09_TPM_calling.r

## call Fold Change, which will be done in R
> Rscript 10_FC.r