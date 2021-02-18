#!/usr/bin/env bash

# copy the extra exercise fastq files and their reference
cp /home/training/chiara/extra_exercise/rCRS.fa ./
cp /home/training/chiara/extra_exercise/fastq/*.fastq.gz ./


# index the fasta of the reference genome
bwa index rCRS.fa


for R1_read in $(ls *R1.fastq.gz);do
	# replace R1 with R2 to get name of paired read
	R2_read="$(echo $R1_read | sed 's/R1/R2/g')" 

	# align the reads
	bwa mem rCRS.fa $R1_read $R2_read | samtools view -b - | samtools sort - -o $R1_read.bam
done
