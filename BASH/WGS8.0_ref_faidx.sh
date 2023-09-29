#!/bin/bash
#SBATCH --ntasks=2 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=faidx 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=faidx.log         
pwd; hostname; date

module load samtools

echo "Indexing FASTA reference genome"

samtools faidx ./certhia.fa

nohup bash WGS8.0_ref_faidx.sh &

date