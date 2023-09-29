#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=sort 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=sort.log         
pwd; hostname; date

module load vcftools

echo "Sorting vcf by position"

vcf-sort ../VC/wgs.wrensc55.vcf > ../VC/wgs.wrensc55.sort.vcf

nohup bash WGS10_sortingvcf.sh &

date