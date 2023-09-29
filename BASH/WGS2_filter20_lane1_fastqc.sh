#!/bin/bash
#SBATCH --ntasks=50 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7gb
#SBATCH --time=4-00:00:00   
#SBATCH --job-name=filter20  
#SBATCH --qos=robinson-b            
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=filter20_%j.log   
pwd; hostname; date

module load fastx_toolkit

echo "Filter sequences with less than 20 Phred scores"

#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/28_S28_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/28_S28_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/28_S28_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/28_S28_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/30_S30_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/30_S30_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/30_S30_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/30_S30_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/17_S17_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/17_S17_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/17_S17_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/17_S17_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/13_S13_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/13_S13_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/13_S13_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/13_S13_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/54_S48_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/54_S48_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/54_S48_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/54_S48_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/37_S36_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/37_S36_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/37_S36_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/37_S36_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/41_S39_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/41_S39_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/41_S39_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/41_S39_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/16_S16_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/16_S16_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/16_S16_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/16_S16_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/29_S29_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/29_S29_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/29_S29_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/29_S29_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/9_S9_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/9_S9_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/9_S9_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/9_S9_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/8_S8_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/8_S8_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/8_S8_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/8_S8_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/20_S20_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/20_S20_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/20_S20_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/20_S20_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/14_S14_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/14_S14_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/14_S14_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/14_S14_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/25_S25_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/25_S25_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/25_S25_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/25_S25_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/19_S19_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/19_S19_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/19_S19_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/19_S19_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/42_S40_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/42_S40_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/42_S40_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/42_S40_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/44_S41_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/44_S41_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/44_S41_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/44_S41_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/27_S27_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/27_S27_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/27_S27_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/27_S27_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/7_S7_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/7_S7_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/7_S7_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/7_S7_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/15_S15_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/15_S15_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/15_S15_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/15_S15_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/47_S42_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/47_S42_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/47_S42_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/47_S42_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/23_S23_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/23_S23_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/23_S23_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/23_S23_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/3_S3_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/3_S3_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/3_S3_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/3_S3_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/57_S50_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/57_S50_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/57_S50_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/57_S50_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/10_S10_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/10_S10_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/10_S10_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/10_S10_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/11_S11_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/11_S11_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/11_S11_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/11_S11_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/1_S1_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/1_S1_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/1_S1_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/1_S1_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/48_S43_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/48_S43_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/48_S43_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/48_S43_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/50_S44_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/50_S44_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/50_S44_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/50_S44_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/18_S18_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/18_S18_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/18_S18_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/18_S18_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/51_S45_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/51_S45_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/51_S45_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/51_S45_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/38_S37_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/38_S37_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/38_S37_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/38_S37_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/40_S38_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/40_S38_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/40_S38_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/40_S38_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/31_S31_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/31_S31_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/31_S31_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/31_S31_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/53_S47_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/53_S47_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/53_S47_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/53_S47_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/12_S12_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/12_S12_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/12_S12_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/12_S12_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/26_S26_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/26_S26_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/26_S26_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/26_S26_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/55_S49_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/55_S49_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/55_S49_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/55_S49_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/4_S4_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/4_S4_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/4_S4_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/4_S4_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/2_S2_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/2_S2_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/2_S2_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/2_S2_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/5_S5_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/5_S5_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/5_S5_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/5_S5_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/24_S24_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/24_S24_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/24_S24_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/24_S24_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/22_S22_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/22_S22_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/22_S22_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/22_S22_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/6_S6_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/6_S6_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/6_S6_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/6_S6_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/35_S34_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/35_S34_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/35_S34_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/35_S34_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/52_S46_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/52_S46_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/52_S46_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/52_S46_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/36_S35_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/36_S35_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/36_S35_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/36_S35_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/32_S32_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/32_S32_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/32_S32_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/32_S32_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/21_S21_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/21_S21_L001_R1_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/21_S21_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/21_S21_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/34_S33_L001_R2_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/34_S33_L001_R2_001.tff.fastq &
fastq_quality_filter -q 20 -p 95 -Q33 -i  /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/34_S33_L001_R1_001.tf.fastq	 -o	/blue/robinson/ldmontalvo/wrens_1/fasta/lane1/34_S33_L001_R1_001.tff.fastq

nohup bash WGS_filter20_lane1_fastqc.sh &

date