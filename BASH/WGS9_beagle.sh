#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=beagle 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=beagle.log         
pwd; hostname; date

module load java

echo "Imputation and QC"

java -jar beagle.r1399.jar gt=all.samples.c5.5.vcf out=wgs.wrensc55.vcf

nohup bash WGS9_beagle.sh &

date