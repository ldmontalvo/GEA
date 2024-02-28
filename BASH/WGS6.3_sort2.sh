#!/bin/bash
#SBATCH --ntasks=5 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb               
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=sort2 
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=sort2.log         
pwd; hostname; date

module load samtools

echo "Sorting by name of the read"



for i in 1	2	#3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	34	35	36	37	38	40	41	42	44	47	48	50	51	52	53	54	55	57	58	59	60	61	62	63	64	65	67	68	69	71	73	75	76	79	80	82	83	85	86	88	91	92	93	94	95	1835	3816	9912	21819	47420	47457	47481	47494	47591	47611	47664	47685	49119	49189	49213	49313	49368	49416	50571	50573	50574	50577	50578
do
# This will sort based on chromosome number and coordinates
samtools sort -o ../Markdup/S${i}.sorted2.bam ../Markdup/S${i}.fix.bam
done

nohup bash WGS6.3_sort2.sh &

date