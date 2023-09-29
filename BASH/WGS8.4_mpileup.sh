#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=mpileup 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=mpileup.log         
pwd; hostname; date

module load samtools
module load bcftools

echo "Variant Calling"


# Variant Calling

samtools mpileup -B -ugf \

# Reference Genome
../Certhia/GCA_018697195.1_ASM1869719v1_genomic.fna \ 

# Combined BAM file
../Markdup/all.samples.merged2.bam | \ 

bcftools call -vmO v -> wgs.wrens.vcf

nohup bash WGS8.2_mpileup.sh &

date