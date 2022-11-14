# EBI Genome Bioinformatics: Scaling Things Up

This is the code repository used for the "Scaling Things Up" section of the EBI course [Genome Bioinformatics](https://www.ebi.ac.uk/training/events/genome-bioinformatics-resequencing-and-variant-calling-2022), named in previous years as "NGS Bioinformatics".

This sections follows the previous 3 days of the course, where command line tools and basic bioinformatics commands to index files and align fastqs to a reference genome have been acquired.
Here we focus on reusing the commands learnt during previous days, to run the same commands using parallelisation and job scheduling.

The following README is a copy of the 2021 Google Docs walkthrough of the interactive part of the session.

## Parallelisation

1. Run git clone on this repository
2. Go into the folder you just cloned, and then inside the “Parallelisation” folder
3. Open the align_all_extra_fqs.sh script. What do you think the script will do?
4. Do you think the script will take a long time to run? What command could we use to time how long a script takes?

![time result for non-parallel alignments](https://user-images.githubusercontent.com/27174192/201706948-f5d7220d-67b2-4be8-a360-d58af19bbe4b.png)

5. Modify the script so that instead of running each alignment, it echos the align command to a file we will call align_commands.sh


*If you try but still can’t do it, you can use the correction in the align_all_extra_fqs__correction.sh file, from the same repository*

6. Run the script using the parallel command, you can even use the time command to measure how long it takes to run

```{bash}
parallel < align_commands.sh
```

7. How long did it take when using parallel to run the command?


## Job Schedulers

1. Remove the `echo` we added to align_all_extra_fqs.sh so that it will run everything in a for loop
2. Do you remember how to submit a job with slurm? (hint: its the `sbatch` command followed by what you want to run)

```{bash}
sbatch align_all_extra_fqs.sh
```

3. Run `squeue` to see your job running. You should see something like this:
![squeue result for non-parallel alignments](https://user-images.githubusercontent.com/27174192/201708540-d36c14be-61b9-4f6c-a366-0771fb656535.png)

4. We will now kill our job, we do this using the `scancel` command followed by the JOBID. For me, this is `scancel 8` . Find your jobid with `squeue` and cancel the job

5. Remove the bam files we generated here

```bash
rm -f *.bam
```

6. Edit the align_all_extra_fqs.sh file to submit each `bwa mem` command to slurm

*This means you wrap the bwa mem line in quotes, and prefix with* sbatch --wrap

```bash
sbatch --wrap "bwa mem rCRS.fa $R1_read $R2_read | samtools view -b - | samtools sort - -o $R1_read.bam"
```

7. See all the jobs running at once
![squeue result for parallel alignments](https://user-images.githubusercontent.com/27174192/201708807-55f1de09-8a0a-48be-8ea6-09b334fa3912.png)
