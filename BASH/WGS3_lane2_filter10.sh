#!/bin/bash
#SBATCH --ntasks=50  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb               
#SBATCH --time=4-00:00:00   
#SBATCH --job-name=lane2_filter10 
#SBATCH --qos=robinson-b              
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=lane2_filter10.log         
pwd; hostname; date

module load parallel fastx_toolkit

echo "Filter sequences with less than 10"


##Eliminate sequences where there is a sinlge Phred score below 10 and then sequences where 5% of reads have a with Phred quality scores below 20
# label of "tff" indicates that the second round of quality filtering is completed
#100% of the bases in a sequence must have a score of higher than 10 for the sequence to be kept
parallel < files_lane1.txt

#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
#fastq_quality_filter -q 20 -p 95 -Q33 -i /ufrc/robinson/ldmontalvo/wrens_1/filter/AS_LDM1_tf.fastq -o /ufrc/robinson/ldmontalvo/wrens_1/filter/AS_LDM1_tff.fastq &


nohup bash WGS3_lane2_filter10.sh &

date