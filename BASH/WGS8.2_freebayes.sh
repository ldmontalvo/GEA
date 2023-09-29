#!/bin/bash
#SBATCH --ntasks=4 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=freebayes
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=freebayes.log         
pwd; hostname; date

module load freebayes/1.3.2

echo "Variant Calling"


# Variant Calling

ls ../Markdup/S*.markdup.rg.bam > bams.list

freebayes -f ../VC/GCA_018697195.1_ASM1869719v1_genomic.fna -m 60 -q 30 -C 5 -! 5 -n 4 -H  --bam-list bams.list > ../VC/all.samples.c5.5.vcf

#-m 60 -q 30 -C 1 -! 1 -n 4 -H
# -m --min-mapping-quality Q Exclude alignments from analysis if they have a mapping quality less than Q.
# -q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q.
# -C --min-alternate-count N Require at least this count of observations supporting. Number of samples that carry an allel 
# an alternate allele within a single individual in order to evaluate the position.
# -! --min-coverage N Require at least this coverage to process a site. How much coverage a sample has to have to be included
# -n --use-best-n-alleles N Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores.
# -H --harmonic-indel-quality Use a weighted sum of base qualities around an indel, scaled by the distance from the indel.

nohup bash WGS8.2_freebayes.sh &

date