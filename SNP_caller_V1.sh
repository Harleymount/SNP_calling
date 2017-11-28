#!/bin/bash

#put reference in FA format with your read 1 and read 2  fastq in a directory that you provide with --directory 
#requirements novoalign, 2bwt-builder, SOAP, Picard, FreeBayes, 
for arg in "$@"; do
  shift
  case "$arg" in
    "--reference") set -- "$@" "-r" ;;
    "--read1") set -- "$@" "-x" ;;
    "--read2")   set -- "$@" "-y" ;;
	"--directory")   set -- "$@" "-z" ;;
	"--samplename")   set -- "$@" "-a" ;;
    *)        set -- "$@" "$arg"
  esac
done


while getopts r:x:y:z:a: option
do
 case "${option}"
 in
 r) reference=${OPTARG};;
 x) r1=${OPTARG};;
 y) r2=${OPTARG};;
 z) dir=$OPTARG;;
 a) sample_name=$OPTARG;;
 esac
done




#---------Process reads-----------------
cd $dir;
#----------trim reads for quality
trimmomatic PE $r1 $r2  output_forward_paired.fastq output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq ILLUMINACLIP:/usr/local/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36;
rm output_forward_unpaired.fastq;
rm output_reverse_unpaired.fastq;

#-----------------map reads to reference using Bowtie 
bowtie2-build $reference my_reference;
bowtie2 --phred33 -p 4 --local -t -x my_reference -X 2000 -1 output_forward_paired.fastq -2 output_reverse_paired.fastq -S $sample_name.sam;
Picard SamFormatConverter I=$sample_name.sam O=$sample_name.bam;
rm *.sam;
Picard AddOrReplaceReadGroups I=$sample_name.bam O=$sample_name-sorted.bam SORT_ORDER=coordinate RGLB=NA RGPL=Illumina RGPU=NA RGSM=$sample_name;
rm $sample_name.bam;
Picard BuildBamIndex I=$sample_name-sorted.bam O=$sample_name-sorted.bai;
Picard MarkDuplicates I=$sample_name-sorted.bam O=$sample_name-final_bowtie.bam M=metrics_bowtie.txt;
Picard BuildBamIndex I=$sample_name-final_bowtie.bam O=$sample_name-final_bowtie.bai;
#=================================================

#-----------------map reads to reference using novoalign 
novoindex reference.nix $reference;
novoalign -d reference.nix -f output_forward_paired.fastq output_reverse_paired.fastq -o SAM > $sample_name.sam 2> log.txt
Picard SamFormatConverter I=$sample_name.sam O=$sample_name.bam;
rm *.sam;
Picard AddOrReplaceReadGroups I=$sample_name.bam O=$sample_name-sorted.bam SORT_ORDER=coordinate RGLB=NA RGPL=Illumina RGPU=NA RGSM=$sample_name;
rm $sample_name.bam;
Picard BuildBamIndex I=$sample_name-sorted.bam O=$sample_name-sorted.bai;
Picard MarkDuplicates I=$sample_name-sorted.bam O=$sample_name-final_novoalign.bam M=metrics_novoalign.txt;
Picard BuildBamIndex I=$sample_name-final_novoalign.bam O=$sample_name-final_novoalign.bai;
#=====================================================================

#-----------------map reads to reference using GSNAP 
gmap_build -d reference $reference;
GSNAP -d reference -A sam output_forward_paired.fastq output_reverse_paired.fastq -t 8 > $sample_name.sam
Picard SamFormatConverter I=$sample_name.sam O=$sample_name.bam;
rm *.sam;
Picard AddOrReplaceReadGroups I=$sample_name.bam O=$sample_name-sorted.bam SORT_ORDER=coordinate RGLB=NA RGPL=Illumina RGPU=NA RGSM=$sample_name;
rm $sample_name.bam;
Picard BuildBamIndex I=$sample_name-sorted.bam O=$sample_name-sorted.bai;
Picard MarkDuplicates I=$sample_name-sorted.bam O=$sample_name-final_GSNAP.bam M=metrics_GSNAP.txt;
Picard BuildBamIndex I=$sample_name-final_GSNAP.bam O=$sample_name-final_GSNAP.bai;
#=====================================================================
rm $sample_name-sorted.bai;
rm $sample_name-sorted.bam;


#==========Now Call Snps for all three alignments============================
#------------------bowtie2 SNP calling
freebayes --fasta-reference $reference $sample_name-final_bowtie.bam --ploidy 1 > $sample_name-bowtie.vcf;
vcffilter -f 'QUAL > 10 & DP > 20  & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' $sample_name-bowtie.vcf > $sample_name-filtered-bowtie2.vcf;
#SNP quality has to be over 30 (QUAL), coverage over 20 (DP)
#also make sure that each read contributes to quality by about 10% (QUAL/AO >10 ) this one is suggested by freebayes presentation im https://www.google.ca/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0ahUKEwilo8ntj9_XAhWc14MKHerhAasQFgguMAE&url=https%3A%2F%2Fwiki.uiowa.edu%2Fdownload%2Fattachments%2F145192256%2Ferik%2520garrison%2520-%2520iowa%2520talk%25202.pdf%3Fapi%3Dv2&usg=AOvVaw0G6VgcVVuS42Bk2WBlP1IS
#and a few other suggested filters like having reads supporting the SNP on both strands (SAF > 0 & SAR > 0)
#and at least two reads balanced on either side of the mutation RPR and RPL > 1


#------------------novoalign SNP calling
freebayes --fasta-reference $reference $sample_name-final_novoalign.bam --ploidy 1 > $sample_name-novoalign.vcf;
vcffilter -f 'QUAL > 10 & DP > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' $sample_name-novoalign.vcf > $sample_name-filtered-novoalign.vcf;

#------------------SOAP SNP calling
freebayes --fasta-reference $reference $sample_name-final_GSNAP.bam --ploidy 1 > $sample_name-GSNAP.vcf;
vcffilter -f 'QUAL > 10 & DP > 20  & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' $sample_name-GSNAP.vcf > $sample_name-filtered-GSNAP.vcf;

mkdir unfiltered_SNPs
mv $sample_name-novoalign.vcf unfiltered_SNPs
mv $sample_name-bowtie.vcf unfiltered_SNPs
mv $sample_name-GSNAP.vcf unfiltered_SNPs

echo 'DONE SNP calling!'
#=================================================

#==============================Now we want to do SNP comparison between the three callers to identify high quality SNPs identified in 2/3 aligners

#-----Add some python/R script that loads all three vcf files and then outputs a single VCF that only contains SNPs shared by 2/3. Then we can do some filtering of these SNPs after



Rscript SNP_filter_V2.R $dir



#a seperate script will be used to compare two SNP files

