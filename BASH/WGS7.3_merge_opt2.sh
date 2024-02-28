#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=Merge 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=Merge.log         
pwd; hostname; date

module load samtools

echo "Merging the BAM files"

# Combining all sorted BAM files
samtools merge ../Markdup/all.samples.merged.bam ../Markdup/*.rmdup.bam

samtools sort all.samples.merged.bam all.samples.merged-sorted.bam

samtools index all.samples.merged-sorted.bam

nohup bash WGS7.2_merge_opt2.sh &

date