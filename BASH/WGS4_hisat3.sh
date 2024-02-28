#!/bin/bash
#SBATCH --ntasks=10 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=Hisat2 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=Hisat2.log         
pwd; hostname; date

module load hisat2
module load parallel

echo "Using Parallel to Assemble Genomes with Hisat2"

parallel --compress < files_hisat2_par.txt

nohup bash WGS4_hisat3.sh &

date