#!/bin/bash
#SBATCH --ntasks=101  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=4-00:00:00   
#SBATCH --job-name=lane2_filter20 
#SBATCH --qos=robinson-b              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=lane2_filter20.log         
pwd; hostname; date

module load parallel fastx_toolkit

echo "Filter sequences with less than 20"


parallel --compress < files_lane2.txt

#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
#fastq_quality_filter -q 20 -p 95 -Q33 -i /ufrc/robinson/ldmontalvo/wrens_1/filter/AS_LDM1_tf.fastq -o /ufrc/robinson/ldmontalvo/wrens_1/filter/AS_LDM1_tff.fastq &


nohup bash WGS3_lane2_filter20.sh &

date