#!/bin/bash
#SBATCH --ntasks=25 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb               
#SBATCH --time=4-00:00:00 
#SBATCH --account=robinson
#SBATCH --qos=robinson-b
#SBATCH --job-name=Hisat2-build
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=Hisat2build.log         
pwd; hostname; date

module load hisat2


cat ./FASTA/chr10.fna \
./FASTA/chr11.fna \
./FASTA/chr12.fna \
./FASTA/chr13.fna \
./FASTA/chr14.fna \
./FASTA/chr15.fna \
./FASTA/chr17.fna \
./FASTA/chr18.fna \
./FASTA/chr1A.fna \
./FASTA/chr1.fna \
./FASTA/chr20.fna \
./FASTA/chr21.fna \
./FASTA/chr23.fna \
./FASTA/chr24.fna \
./FASTA/chr2.fna \
./FASTA/chr3.fna \
./FASTA/chr4A.fna \
./FASTA/chr4.fna \
./FASTA/chr5.fna \
./FASTA/chr6.fna \
./FASTA/chr7.fna \
./FASTA/chr8.fna \
./FASTA/chr9.fna \
./FASTA/chrZ.fna \
./FASTA/chr19.fna > certhia.fa



#hisat2-build -f -p 25 certhia.fa ./Indexes/Certhia


nohup bash WGS4_hisat2-build.sh &

date

