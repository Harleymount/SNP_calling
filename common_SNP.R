#this script will be used to compare the SNPs generated using the SNP_Caller_V1.sh script
#input format name should be high_quality_snps_samplename.vcf where you adjust the samplename part for each 

library(stringi)
args = commandArgs(trailingOnly=TRUE)
vcf1=args[1]
vcf2=args[2]

vcf1<-'/Users/Harley/Desktop/SNP/high_quality_SNPS_1.csv'
vcf2<-'/Users/Harley/Desktop/SNP/high_quality_SNPS_2.csv'

sample_name_1<-stri_sub(unlist(strsplit(vcf1,'_'))[4],1,-5)
sample_name_2<-stri_sub(unlist(strsplit(vcf2,'_'))[4],1,-5)

sample1<-read.csv(vcf1)
sample2<-read.csv(vcf2)


shared_snps_sample_1<-sample1[sample1$unique_SNP_ID %in% sample2$unique_SNP_ID,]
shared_snps_sample_2<-sample2[sample2$unique_SNP_ID %in% sample1$unique_SNP_ID,]
shared_snps_sample_1<-cbind(shared_snps_sample_1, sample=sample_name_1)
shared_snps_sample_2<-cbind(shared_snps_sample_2, sample=sample_name_2)

merged<-data.frame()

merged<-rbind(shared_snps_sample_1,shared_snps_sample_2)
merged$X<-NULL
col_names<-c('contig','position','ref','alternate','qual','info','format','ploidy_info','contig', 'mapper','unique_SNP_ID', 'sample_name')
colnames(merged)<-col_names 
write.csv(merged, 'common_SNPs.csv')
