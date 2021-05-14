# RNA-seq_data_processing

## load the packages
> module load SRA-Toolkit FastQC Trimmomatic hisat2 SAMtools GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA HTSeq
#### Note that if the fastq-dump doesn't work, run the code below to save the config for sra toolkit
> vdb-config --interactive  

## build the genome index
> hisat2-build Pvirgatum_516_v5.0.fa Pvirgatum_516_v5.0.fa_gi

## You can apply the code 01_RNA-seq_processing.py to do the jobs automatically. An example command line for this code is as follow, which will produce a .sh file to be submitted to the slurm queue

> python 01_RNA_seq_processing.py -SRA SRR14066076 -genome_seq Pvirgatum_516_v5.0.fa -gff Pvirgatum_516_v5.1.gene.gff3 -layout PE -workdir /mnt/scratch/peipeiw/Data_for_Kenia/For_RNA_seq_pipeline -trim y -adapters all_PE_adapters.fa

## combine the read count files
> python 03_combine_read_counts.py

## call the TPM, which will be done in R
> Rscript 04_TPM_calling.r Read_counts.txt Pvirgatum_516_v5.1_transcript_length.txt
