#!/bin/bash
#SBATCH --ntasks=12 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb               
#SBATCH --time=4-00:00:00
#SBATCH --account=robinson
#SBATCH --qos=robinson-b
#SBATCH --job-name=Hisat2 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=Hisat2.log         
pwd; hostname; date

module load hisat2

echo "Indexing the reference genome"


hisat2-build -f reference.genome.fna Certhia


echo "Assembling to Reference Genome"



for i in 1	2 
do
hisat2 -p 12 -q --phred33 \
--rg-id S${i} --rg SM:S${i} --rg PL:illumina \
-x ./Indexes/Certhia \
-1 ./${i}_R1_001.tff.fastq.gz \
-2 ./${i}_R2_001.tff.fastq.gz \
-S ./S${i}.sam
done

nohup bash WGS4_hisat2.sh &

date

