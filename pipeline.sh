#!/usr/bin/env bash
PATH_TO_REF_GENOME=/data/extra_exercise_EBI_NGS_sept2019/rCRS.fa
PATH_TO_READS=/data/extra_exercise_EBI_NGS_sept2019/fastq/
NB_THREADS=4

# set variables to our programs not in $PATH
trimmomatic="java -jar /home/software/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
picard="java -jar /home/software/picard/picard.jar"
gatk="/home/software/gatk/gatk-4.1.3.0/gatk"


# Run trimmomatic on the input fastq files
output_dir=/data/output_files/1_trimmomatic_output

if [ ! -d "$output_dir" ]; then
	mkdir -p $output_dir

	read_ids=`ls $PATH_TO_READS | sed 's/\(_R\).*/\1/g' | sort -u`
	for read_id in $read_ids
	do
		echo $read_id
		read_1=`paste <(echo "$read_id") <(echo "1.fastq.gz") --delimiters ''`
		read_2=`paste <(echo "$read_id") <(echo "2.fastq.gz") --delimiters ''`
		$trimmomatic PE -phred33 -threads $NB_THREADS -trimlog $output_dir/$read_id \
			$PATH_TO_READS/$read_1 $PATH_TO_READS/$read_2 \
			$output_dir/$read_1 $output_dir/$read_1.unpaired \
			$output_dir/$read_2 $output_dir/$read_2.unpaired \
			LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15  MINLEN:36
	done
fi


# samtools index ref genome
if [ ! -f "$PATH_TO_REF_GENOME.fai" ]; then
	cd $(dirname $PATH_TO_REF_GENOME) && samtools faidx $PATH_TO_REF_GENOME
fi

# bwa index ref genome
if [ ! -f "$PATH_TO_REF_GENOME.bwt" ]; then
	cd $(dirname $PATH_TO_REF_GENOME) && bwa index -a is $PATH_TO_REF_GENOME
fi

# create seq dictionary
if [ ! -f "$(dirname $PATH_TO_REF_GENOME)/picard_dict.dict" ]; then
	cd $(dirname $PATH_TO_REF_GENOME) && $picard CreateSequenceDictionary R=$PATH_TO_REF_GENOME O=picard_dict.dict
fi


# Use bwa mem for mapping
output_dir=/data/output_files/2_mapping_output

if [ ! -d "$output_dir" ]; then
	mkdir -p $output_dir
	read_ids=`ls $PATH_TO_READS | sed 's/\(_R\).*/\1/g' | sort -u`
	for read_id in $read_ids
	do
		read_1=`paste <(echo "$read_id") <(echo "1.fastq.gz") --delimiters ''`
		read_2=`paste <(echo "$read_id") <(echo "2.fastq.gz") --delimiters ''`

		bwa mem -t $NB_THREADS -R '@RG\tID:'"$n"'\tLB:library\tPL:Illumina\tPU:lane'"$n"'\tSM:yeast' \
		$PATH_TO_REF_GENOME \
		$PATH_TO_READS/$read_1 $PATH_TO_READS/$read_2 > $output_dir/$read_id.unsorted.bam

		samtools sort -@ $NB_THREADS $output_dir/$read_id.unsorted.bam $output_dir/$read_id.sorted.bam
		cd $output_dir && samtools index $read_id.sorted.bam && cd ..

		#$picard MarkDuplicates INPUT=$output_dir/$read_id.sorted.bam OUTPUT=$output_dir/$read_id.deduped.sorted.bam METRICS_FILE=$output/$read_id.dupl_metrics.txt
	done
fi


# Use GATK for variant calling







