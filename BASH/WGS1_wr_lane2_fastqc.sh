#!/bin/bash
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=30-00:00:00   
#SBATCH --job-name=check_test2             
#SBATCH --mail-type=END,FAIL               
#SBATCH --mail-user=ldmontalvo@ufl.edu                    
#SBATCH --output=check_test2_%j.log         
pwd; hostname; date

module load fastqc

echo "check quality of data using FastQC"


fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_82_83_S778_S522_L002/83_S19_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_82_83_S778_S522_L002/83_S19_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_106_50577_S720_S514_L002/50577_S38_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_106_50577_S720_S514_L002/50577_S38_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_74_75_S721_S508_L002/75_S14_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_74_75_S721_S508_L002/75_S14_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_100_47481_S744_S552_L002/47481_S32_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_100_47481_S744_S552_L002/47481_S32_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_105_47611_S741_S584_L002/47611_S37_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_105_47611_S741_S584_L002/47611_S37_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_75_76_S723_S502_L002/76_S15_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_75_76_S723_S502_L002/76_S15_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_97_47664_S765_S549_L002/47664_S29_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_97_47664_S765_S549_L002/47664_S29_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_66_67_S767_S572_L002/67_S9_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_66_67_S767_S572_L002/67_S9_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_58_59_S775_S509_L002/59_S2_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_58_59_S775_S509_L002/59_S2_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_87_88_S792_S534_L002/88_S22_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_87_88_S792_S534_L002/88_S22_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_99_47420_S719_S569_L002/47420_S31_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_99_47420_S719_S569_L002/47420_S31_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_109_47685_S787_S588_L002/47685_S41_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_109_47685_S787_S588_L002/47685_S41_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_91_92_S701_S597_L002/92_S24_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_91_92_S701_S597_L002/92_S24_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_98_49368_S793_S573_L002/49368_S30_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_98_49368_S793_S573_L002/49368_S30_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_92_93_S754_S535_L002/93_S25_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_92_93_S754_S535_L002/93_S25_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_61_62_S724_S547_L002/62_S5_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_61_62_S724_S547_L002/62_S5_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_115_47494_S733_S561_L002/47494_S46_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_115_47494_S733_S561_L002/47494_S46_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_116_9912_S762_S512_L002/9912_S47_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_116_9912_S762_S512_L002/9912_S47_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_110_49119_S703_S539_L002/49119_S42_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_110_49119_S703_S539_L002/49119_S42_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_85_86_S740_S570_L002/86_S21_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_85_86_S740_S570_L002/86_S21_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_103_49416_S748_S553_L002/49416_S35_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_103_49416_S748_S553_L002/49416_S35_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_79_80_S745_S557_L002/80_S17_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_79_80_S745_S557_L002/80_S17_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_81_82_S746_S583_L002/82_S18_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_81_82_S746_S583_L002/82_S18_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_67_68_S738_S525_L002/68_S10_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_67_68_S738_S525_L002/68_S10_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_118_49189_S713_S586_L002/49189_S48_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_118_49189_S713_S586_L002/49189_S48_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_108_49313_S789_S579_L002/49313_S40_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_108_49313_S789_S579_L002/49313_S40_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_104_50574_S768_S558_L002/50574_S36_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_104_50574_S768_S558_L002/50574_S36_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_57_58_S743_S536_L002/58_S1_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_57_58_S743_S536_L002/58_S1_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_111_47591_S776_S537_L002/47591_S43_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_111_47591_S776_S537_L002/47591_S43_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_113_50573_S730_S580_L002/50573_S45_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_113_50573_S730_S580_L002/50573_S45_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_84_85_S781_S510_L002/85_S20_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_84_85_S781_S510_L002/85_S20_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_112_50571_S782_S528_L002/50571_S44_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_112_50571_S782_S528_L002/50571_S44_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_78_79_S707_S596_L002/79_S16_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_78_79_S707_S596_L002/79_S16_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_102_50578_S718_S515_L002/50578_S34_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_102_50578_S718_S515_L002/50578_S34_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_95_1835_S717_S595_L002/1835_S28_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_95_1835_S717_S595_L002/1835_S28_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_93_94_S798_S527_L002/94_S26_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_93_94_S798_S527_L002/94_S26_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_90_91_S705_S563_L002/91_S23_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_90_91_S705_S563_L002/91_S23_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_70_71_S711_S560_L002/71_S12_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_70_71_S711_S560_L002/71_S12_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_64_65_S726_S577_L002/65_S8_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_64_65_S726_S577_L002/65_S8_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_94_95_S751_S542_L002/95_S27_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_94_95_S751_S542_L002/95_S27_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_107_3816_S731_S571_L002/3816_S39_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_107_3816_S731_S571_L002/3816_S39_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_62_63_S790_S594_L002/63_S6_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_62_63_S790_S594_L002/63_S6_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_68_69_S766_S555_L002/69_S11_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_68_69_S766_S555_L002/69_S11_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_63_64_S791_S529_L002/64_S7_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_63_64_S791_S529_L002/64_S7_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_60_61_S706_S574_L002/61_S4_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_60_61_S706_S574_L002/61_S4_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_120_49213_S709_S575_L002/49213_S50_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_120_49213_S709_S575_L002/49213_S50_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_72_73_S759_S516_L002/73_S13_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_72_73_S759_S516_L002/73_S13_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_59_60_S750_S532_L002/60_S3_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_59_60_S750_S532_L002/60_S3_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_101_47457_S785_S504_L002/47457_S33_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_101_47457_S785_S504_L002/47457_S33_L002_R1_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_119_21819_S736_S543_L002/21819_S49_L002_R2_001.fastq.gz
fastqc /blue/robinson/ldmontalvo/wrens_1/NS2648-SRobinson_S4-HHGHGDSX3-Lane2/NS2648_119_21819_S736_S543_L002/21819_S49_L002_R1_001.fastq.gzd


date