#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=filter_vcf 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=filter_vcf.log         
pwd; hostname; date

module load vcftools
#module load plink

echo "Filtering vcf"

vcftools --vcf ../VC/wgs.wrensc55.sort.vcf --max-missing 0.8 --maf 0.05 --hwe 0.01 --recode --recode-INFO-all --out ../ImputationQC/wgs.wrensc55.filter.vcf
#--hap-r2 --min-r2 .7 --ld-window-bp 500000 \
#--out ../VC/LDs

echo "Making PLINK files"

#plink --vcf ../VC/wgs.wrens.filter.vcf --make-bed --out ../ImputationQC/wgs.wrens.filter

nohup bash WGS11_filter_vcf.sh &

date