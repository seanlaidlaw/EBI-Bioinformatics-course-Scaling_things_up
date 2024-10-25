#!/usr/bin/env bash
set -e # stop script on error

# copy the extra exercise fastq files and their reference
cp /home/training/scaling_things_up/chr3.* ./
cp /home/training/scaling_things_up/*.fastq.gz ./


for R1_read in $(ls *R1.fastq.gz);do
	# replace R1 with R2 to get name of paired read
	R2_read="$(echo $R1_read | sed 's/R1/R2/g')" 

	# align the reads
	bwa mem chr3.fa $R1_read $R2_read | samtools view -b - | samtools sort - -o $R1_read.bam
done