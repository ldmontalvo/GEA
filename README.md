# NOTES ON WGS ANALYSES FOR WRENS
## By Luis Daniel Montalvo

This repository contains data quality control, read mapping, and genome-environment association analyses for 100 Campylorhynchus individuals sequenced with whole-genome resequencing. The sampled individuals span an environmental gradient of climate variables. The workflows and scripts provided here allow for quality filtering of raw reads, alignment to the reference genome, identification of climate-associated SNPs, and functional characterization of candidate adaptive genes. 

## Downloading the Reference Genome

1.	We will use the genome of Certhia americana (ASM1869719v1) as the reference genome (https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_018697195.1/). If we click on “FTP directory for GenBank assembly” we will see the different files we can download including the sequences and annotations.
2.	We create a directory where we can download the files. In this case, my folder is “Reference”.
3.	We use the following command line to download the files into our folders.
rsync --copy-links --times --verbose --recursive rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/697/195/GCA_018697195.1_ASM1869719v1/ ./Reference
More information about downloading NCBI genomes: https://www.youtube.com/watch?v=-X0H024ST8k
4.	Unzip (decompress) the Fasta file that could come with the extension “.fna”.
gunzip GCA_018697195.1_ASM1869719v1_genomic.fna.gz
5.	Creating the indexes for the reference genome. The process took ~10 minutes:
module load hisat2
hisat2-build -p 10 -f ./GCA_018697195.1_ASM1869719v1_genomic.fna Certhia_hisat2
6.	It could be useful to change the file names to make easier the use of loop in bash commands. The lines below could be used for that purpose.
for i in *; do mv $i ${i/S*002_/}; done
this code remove everything (*) from “S” to “002_” and replace it with nothing  (whatever between last “/” and “}”).

7.	I deleted sequences where there is a single Phred score below 10. Then I filtered sequences where 5% of reads have a Phred quality scores below 20 (using fastq_quality_filter).
module load fastx_toolkit
echo "Filter sequences with less than 10 and 20 Phred scores"
#95% of the bases in a sequence must have a score of more than 20 for the sequence to be kept
fastq_quality_filter -q 20 -p 95 -Q33 -i /blue/robinson/ldmontalvo/wrens_1/fasta/lane1/*.tf.fastq &
nohup WGS1_wr_lane1_fastqc.sh &
date

## Genome Assembly
8.	I used hisat2. I had the issue that the software took the RefId instead of the chr name. I asked Sarah and she also had this issue. She said some people just continue with the notation of RefID instead of chr name. I could change these names later. I used the script called WGS4_hisat2. 
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

9.	I converted Sam files from hisat2 into Bam using the WGS5_bam script. 
module load samtools
echo "Converting to bam"
for i in 1	2	# other samples here
do
samtools view -S -b -@ 4 ./S${i}.sam -o ./S${i}.bam 
samtools index ./S${i}.bam
done
nohup bash WGS5_bam.sh &
date

## Removing PCR Duplicates
10.	I used the function markdup of Samtools. Rmdup was deprecated. I followed the steps describe at https://www.biostars.org/p/415831/. 
11.	First, I had to sort the files by the names of the reads. I used the WGS6_sort.sh script with the following line with a loop.

samtools sort -n -o ../Markdup/S${i}sorted.bam -O BAM ../Bam/S${i}.bam

12.	Then I fill in mate coordinates and insert size fields with the WGS6_fixmate.sh script with the following line with a loop.
samtools fixmate -m ../Markdup/S${i}sorted.bam ../Markdup/S${i}fix.bam
13.	I sort based on chromosome number and coordinates using the WGS6_sort2.sh with the following line with a loop.
samtools sort -o ../Markdup/S${i}sorted.bam ../Markdup/S${i}fix.bam
14.	I removed the duplicates and print some basic stats about the resulting files with WGS6_markdup.sh script with the following line with a loop.
samtools markdup -r -s ../Markdup/S${i}sorted.bam  ../Markdup/S${i}.rmdup.bam
15.	I generated indexes for each BAM file with the following script WGS7.1_sort_indexing_merge.sh.
module load samtools
echo "Indexing and Merging"
for i in 1	2	… # Put the rest of samples here
do
samtools index ../Markdup/S${i}.markdup.bam
done
#samtools merge -r ./all.samples.merged2.bam ./*.markdup.bam
nohup bash WGS7.1_sort_indexing_merge.sh &
date

The script has one line with the command samtools merge. I could not make the -r flags to work. This is to generate a Read Group tag with the sample name at the end of each file. Instead, I used the alternative below.
16.	I wrote RG tags in each BAM files with the following script WGS7.2_rg_tags.sh.
for i in 1	2… # Put the rest of samples here
do
samtools addreplacerg -r "@RG\tID:${i}\tSM:${i}" -o ../Markdup/S${i}.markdup.rg.bam ../Markdup/S${i}.rmdup.bam
#samtools addreplacerg -r "@RG\tID:$name\tPG:samtools addreplacerg\tSM:file_name}"
done
nohup bash rg_tag.sh &
date

17.	We also make a fasta indexed of the reference genome to split it into regions with the following script WGS8.0_ref_faidx.sh.
module load samtools
echo "Indexing FASTA reference genome"
samtools faidx ../Certhia/GCA_018697195.1_ASM1869719v1_genomic.fna
nohup bash WGS9_ref_faidx.sh &
date
Now the BAM files should be ready for the Variant Calling.

## Variant Calling (VC)
### Option 1
18.	We used Freebayes (Garrison and Marth 2012) for variant calling running the script 8.2_freebayes.sh. We used a weighted sum of base qualities around an indel, minimum mapping quality of 60, base quality of 30, minimum coverage of 5x, minimum alternate allele count of 5, and a maximum of 4 best alleles considered.
module load freebayes/1.3.2
echo "Variant Calling"
### Variant Calling
ls ../Markdup/S*.markdup.rg.bam > bams.list
freebayes -f ../VC/GCA_018697195.1_ASM1869719v1_genomic.fna -m 60 -q 30 -C 1 -! 1 -n 4 -H  --bam-list bams.list > ../VC/all.samples.c5.5.vcf
#-m 60 -q 30 -C 1 -! 1 -n 4 -H
#-m --min-mapping-quality Q Exclude alignments from analysis if they have a mapping quality less than Q.
#-q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q.
#-C --min-alternate-count N Require at least this count of observations supporting 
#an alternate allele within a single individual in order to evaluate the position.
#-! --min-coverage N Require at least this coverage to process a site.
#-n --use-best-n-alleles N Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores.
#-H --harmonic-indel-quality Use a weighted sum of base qualities around an indel, scaled by the distance from the indel.
nohup bash WGS8.2_freebayes.sh &
date

### Option 2
19.	For this options we could use Samtools and bcftools with the script WGS8.4_mpileup.sh.
module load samtools
module load bcftools
echo "Variant Calling"
### Variant Calling
samtools mpileup -B -ugf \
### Reference Genome
../Certhia/GCA_018697195.1_ASM1869719v1_genomic.fna \ 
### Combined BAM file
../Markdup/all.samples.merged.bam | \ 
bcftools call -vmO v -> wgs.wrens.vcf
nohup bash WGS9_mpileup.sh &
date

More information at:
Freebayes
https://github.com/freebayes/freebayes
Samtools
https://wikis.utexas.edu/display/bioiteam/Variant+calling+using+SAMtools
https://samtools.github.io/bcftools/howtos/variant-calling.html

## Imputation
20.	Utilizing genotype correction and imputation techniques that make use of all sample data and linkage disequilibrium between variants can significantly increase genotyping accuracy in the case of inadequate sampling or coverage (Nielsen et al. 2011). To genotype samples at variants with missing data and correct for missing genotypes, we used Beagle v.4 (Browning and Browning 2007) in our imputation strategy. This method has been shown to dramatically increase genotyping accuracy with low coverage data (Nielsen et al. 2011). We used the script WGS9_beagle.sh.
module load java
echo "Imputation and QC"
java -jar beagle.r1399.jar gt=all.samples.c5.5.vcf out=wgs.wrensc55.vcf
nohup bash WGS9_beagle.sh &
date

21.	I could make the most recent version of Beagle work (beagle.22Jul22.46e.jar), but the beagle.r1399 worked with any issue.
22.	After imputation, I decompressed and sorted the resulting vcf.gz files. I used vcftools and the script WGS10_sortingvcf.sh. The resulting file of the sorting without decompression first produce a gz file that I couldn’t open in other tools. Decompression before sorting didn’t have that issue.
module load vcftools
echo "Sorting vcf by position"
vcf-sort ../VC/wgs.wrens.vcf.gz > ../VC/wgs.wrensc55.sort.vcf.gz
nohup bash WGS10_sortingvcf.sh &
date

## Quality Control
23.	We filtered the resulting vcf files after imputation with Beagle and sorting. We calculated the missing percentage of data per sample basis with the following command line. This line produces a data.ed.imiss file with the last column having the percentage of missing data for the sample in that row. 

vcftools --vcf ../VC/wgs.wrensc55.sort.vcf --missing-indv --out datac55.ed

24.	The imiss file showed all samples with 0% of missing data. However, I used the flag –max-missing 0.8 in vcftools.  I filtered with a minimum allele frequency (maf) of 0.05, produced a file with the alleles in disequilibrium, and removed alleles not in HWE using a p-value of 0.01. I used the script WGS11_filter_vcf.sh. Marker quality and other filters were already used in the variant calling. This script can estimate the LD and produce a report and filter by HWE, max-missing, and maf if these lines are activated (remove the #).

module load vcftools
echo "Filtering vcf"
vcftools --vcf ../VC/wgs.wrens.sort.vcf \
--hap-r2 --min-r2 .7 --ld-window-bp 500000 \ # Linkage Disequilibrium
--out LDs
#--max-missing 0.8 --maf 0.05 \
#--hwe 0.01 \ # HWE
#--recode --recode-INFO-all --out ../VC/wgs.wrens.filter.vcf 
nohup bash WGS11_filter_vcf.sh &
date

25.	Transform the filtered vcf to PED format. In this conversion, I am using the file GCA_018697195.1_ASM1869719v1_assembly_report.txt in the Certhia reference genome to make a Txt file with the access numbers and the corresponding chromosome name. The file chr.names.txt has the access number in the first column and the chromosome names in the second. With chr.names.txt and the filtered vcf from the previous step, I ran the following line from WGS11_filter_vcf2.sh script.
vcftools --vcf ../ImputationQC/wgs.wrens.filter.vcf --plink --chrom-map ../ImputationQC/chr.names.txt --out ../ImputationQC/wgs.wrensc55.upname

bcftools annotate --rename-chrs ./chr.names.txt wgs.wrensc55.filter.vcf.recode.vcf > wgs.wrensc55.renamedChr.vcf

26.	We got the wgs.wrensc55.upname.ped and wgs.wrensc55.upname.map files from the previous step. We filter now by Linkage Disequilibrium with PLINK. First, we get the estimates of R2 and the list of variants to prune.
plink --ped ../ImputationQC/wgs.wrens.upname.ped --map ../ImputationQC/wgs.wrens.upname.map --aec --indep-pairwise 50 5 0.8 --make-bed --out ../ImputationQC/wgs.wren.ld

27.	We got the wgs.wren.ld.bed, wgs.wren.ld.prune.in and wgs.wren.ld.prune.out files from the previous step.  We used these files to filter LD with the following line in script WGS12_filter_LD.sh.
plink --bfile ../ImputationQC/wgs.wren.ld.bed --aec --exclude ../ImputationQC/wgs.wren.ld.prune.out --allow-extra-chr --recode ped  –out ../ImputationQC/wgs.wrensc55.filteredLD


28.	The command line above should produce already a PLINK flat (ped and map). But I can convert from PLINK binary to PLINK flat with the command line below. TASSEL asks for the flat format.
plink --bfile wgs.wrensc55.filteredLD --allow-extra-chr --recode vcf --out wgs.wrensc55.prunedLD

## Analysis using LEA
29.	Getting a testing subsample, extractin chromosome seven.

grep -E '^(#|2[[:space:]])' wgs.wrens.named.vcf > wgs.wrens.chr7.vcf

30.	The rest of the analyses with LEA are described in the script WGS_LEA.Rmd.
31.	We used these commands to get missing data.

bcftools stats wgs.wrens.named.vcf

vcftools --vcf wgs.wrens.named.vcf  --missing-indv

vcftools --vcf wgs.wrens.named.vcf  --missing-site
I got some files for missing data per sample (out.imiss) and per site (out.lmiss)
32.	We obtained the heterozygosity per sample with the following line.

plink --vcf wgs.wrens.named.vcf --sample-counts cols=hom,het --allow-extra-chr
Blasting Candidates Genes
1.	First I tried the candidates gene for the Radseq data. I’ll work in Visual Studio Code (VSC). I installed the necessary package in Ubuntu terminal in VSC such as Anaconda. Install Anaconda according to the following website.
https://docs.anaconda.com/free/anaconda/install/linux/
2.	I install Blast+ following the steps at:
https://www.ncbi.nlm.nih.gov/books/NBK279690/
3.	I can get a simple search by running:

blastn -query sig.snps.fa -db nt -evalue 1e-20 -out blastn.txt -entrez_query "Taeniopygia guttata [ORGN]" -parse_deflines -remote

4.	I can also run a local query if I download the reference data. Download the Zebra Finch annotated genome in the link below. This is the project GCF_003957565.2 for Taeniopygia guttate (zebra finch) and the genome assembly bTaeGut1.4.pri which is a revised version 4 assembly of the male zebra finch bTaeGut1
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003957565.2/
5.	Unzip the zip file.

unzip ncbi_dataset.zip

6.	Prepare the database before running Blast+.

makeblastdb -in GCF_003957565.2_bTaeGut1.4.pri_genomic.fna -out db2/taegut_prot -dbtype prot 

7.	The file I am working with is GCA_003957565.4_bTaeGut1.4.pri_genomic.fna. Run the following command in the terminal. 

blastn -query sig.snps.fa -db db2/taegut_prot -evalue 1e-20 -out blastn.txt
