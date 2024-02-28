#!/bin/bash
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=check_test1              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test1_%j.log         
pwd; hostname; date

module load fastx_toolkit

echo "Filter sequences with less than 10 and 20 Phred scores"


#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/*.tf.fastq &


nohup bash WGS


date