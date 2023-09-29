#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=ld_plink2 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=ld_plink2.log         
pwd; hostname; date


module load  plink/2.00a3LM

echo "Filtering LD"

#plink --ped ../ImputationQC/wgs.wrensc55.upname.ped --map ../ImputationQC/wgs.wrensc55.upname.map --aec --indep-pairwise 50 5 0.8 --make-bed --out ../ImputationQC/wgs.wrenc55.ld

#--update-name ../ImputationQC/chr.names.txt --set-all-var-ids

plink --bfile ../ImputationQC/wgs.wrenc55.ld --exclude ../ImputationQC/wgs.wrenc55.ld.prune.out --allow-extra-chr  --recode ped --out ../ImputationQC/wgs.wrensc55.filteredLD 

nohup bash WGS12_filter_LD.sh &

date