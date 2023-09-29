#!/bin/bash
#SBATCH --ntasks=16 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=2freebayes
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=2freebayes.log         
pwd; hostname; date

module load freebayes/1.3.2

echo "Variant Calling"


# Variant Calling

ls ./S*.markdup.rg.bam > bams.list1

freebayes-parallel <(./split_ref_by_bai_datasize.py ./GCA_018697195.1_ASM1869719v1_genomic.fna.fai 100000) 15 -f ./GCA_018697195.1_ASM1869719v1_genomic.fna -m 60 -q 30 -C 1 -! 1 -n 4 -H --bam-list bams.list1 > ./all.samples.freeBayes2.vcf

#-m --min-mapping-quality Q Exclude alignments from analysis if they have a mapping quality less than Q.
# -q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q.
# -C --min-alternate-count N Require at least this count of observations supporting 
# an alternate allele within a single individual in order to evaluate the position.
# -! --min-coverage N Require at least this coverage to process a site.
# -n --use-best-n-alleles N Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores.
# -H --harmonic-indel-quality Use a weighted sum of base qualities around an indel, scaled by the distance from the indel.

nohup bash WGS8.3_freebayes_parallel.sh &

date