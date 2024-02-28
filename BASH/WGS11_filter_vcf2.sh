#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=vcf_plink 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=vcf_plink.log         
pwd; hostname; date

module load vcftools
#module load plink

#echo "Filtering vcf"

#vcftools --vcf ../VC/wgs.wrens.sort.vcf --max-missing 0.8 --maf 0.05 --hwe 0.01 --recode --recode-INFO-all --out ../ImputationQC/wgs.wrens.filter.vcf
#--hap-r2 --min-r2 .7 --ld-window-bp 500000 \
#--out ../VC/LDs


echo "Making PLINK files"

#plink --bfile ../ImputationQC/wgs.wrens.filter  --update-map chr.names.txt --update-chr --make-bed --out ../ImputationQC/wgs.wrens.upname
vcftools --vcf ../ImputationQC/wgs.wrensc55.filter.vcf.recode.vcf --plink --chrom-map ../ImputationQC/chr.names.txt --out ../ImputationQC/wgs.wrensc55.upname

nohup bash WGS11_filter_vcf2.sh &

date