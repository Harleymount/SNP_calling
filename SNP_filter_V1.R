#====================Take user input VCFs
args = commandArgs(trailingOnly=TRUE)
dir=args[1]




#-------------------variables for filtering SNPs 
distance_to_end_of_contig=250#in the future take a config file like julio, but for now we can either hard code it in or just change it once here at the top
files <- list.files(path=dir, pattern="*.vcf", full.names=T, recursive=FALSE)
setwd(dir)
library(stringi)




#-------------------get sample names based on VCF files, could be useful in teh future but not really necessary at the moment
sample_names<-c()
for (i in files){
    sample<-stri_sub(tail(unlist(strsplit(i, '/')),n=1), 1, -5)
    sample_names<-c(sample_names, sample)
}




#--------------------SNP filtering in R where VCFfilter did not filter
#for now I will hardcode it in so that it always checks for bowtie, GSNAP and novoalign so I will harcode in dataframes for these samples and populate them
#read in VCFs and filter those near the ends of contigs




#--------------------Bowtie2 
bowtie2_df<-read.table(files[1])
bowtie2_filtered_df<-data.frame()
for (i in 1:length(bowtie2_df[,1])){
    length<-as.integer(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[4]) 
    coverage<-as.integer(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[6])
    identifier<-cat(as.character(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[1]),'_',as.character(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[2]),'_',as.character(bowtie2_df$V2[i]))
    position<-as.integer(bowtie2_df$V2[i])
    if(position  > distance_to_end_of_contig & (length-position > distance_to_end_of_contig)){
        bowtie2_filtered_df<-rbind(bowtie2_filtered_df, bowtie2_df[i,])
    }
    
}
#--------------------------------------------------------------
#------------------------------GSNAP
GSNAP_df<-read.table(files[2])
GSNAP_filtered_df<-data.frame()
for (i in 1:length(GSNAP_df[,1])){
    length<-as.integer(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[4]) 
    coverage<-as.integer(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[6])
    identifier<-cat(as.character(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[1]),'_',as.character(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[2]),'_',as.character(GSNAP_df$V2[i]))
    position<-as.integer(GSNAP_df$V2[i])
    if(position  > distance_to_end_of_contig & (length-position > distance_to_end_of_contig)){
        GSNAP_filtered_df<-rbind(GSNAP_filtered_df, GSNAP_df[i,])
    }
    
}
#-----------------------------------------------------------------
#--------------------------------novoalign
novoalign_df<-read.table(files[3])
novoalign_filtered_df<-data.frame()
for (i in 1:length(novoalign_df[,1])){
    length<-as.integer(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[4]) 
    coverage<-as.integer(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[6])
    identifier<-cat(as.character(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[1]),'_',as.character(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[2]),'_',as.character(novoalign_df$V2[i]))
    position<-as.integer(novoalign_df$V2[i])
    if(position  > distance_to_end_of_contig & (length-position > distance_to_end_of_contig)){
        novoalign_filtered_df<-rbind(novoalign_filtered_df, novoalign_df[i,])
    }
    
}
#-------------------------------------------------------------------

#-------------Now want to filter SNPs found in a cluster, this will be loosely defined as more than 3 SNPs within a 10 bp window 
#the filtration for clusters comes from here : Ren, C. et al. SNP discovery by high-throughput sequencing in soybean. BMC Genomics 11, 469 (2010).

#get the identifier --> node and position for all SNPs and then filter within tables first to remove stretches of SNPs 
#then split SNPs into high quality (found in all three) or  medium qualtiy (found in 2/3) or low quality (found in 1/3) samples





#---------------------------------------------SNP cluster program for finding clusters of SNPs and excluding them (likely mapping errors)

#------------Bowtie 
#get bowtie SNPs
contig_list<-c()
for (i in 1:length(bowtie2_filtered_df[,1])){
    contig_info<-unlist(strsplit(as.character(bowtie2_filtered_df$V1[i]), '_'))
    contig_identifier<-paste0(contig_info[1],'_',  contig_info[2])
    contig_list<-c(contig_list, contig_identifier)
}


#add the identifier column onto the dataframe so that we can use these names to filter out later 
bowtie2_filtered_df<-cbind(bowtie2_filtered_df, contig_list)
#now that we have the SNPs with unique identifiers we can compare every SNP in the table to the others to see if there are clusters

SNP_clusters<-c()
for (i in 1:length(bowtie2_filtered_df[,1])){
    position<-as.integer(bowtie2_filtered_df$V2[i])
    node<-bowtie2_filtered_df$contig_list[i]
    query_table<-subset(bowtie2_filtered_df, bowtie2_filtered_df$contig_list==node)
    upper_bound=as.integer(position)+5
    lower_bound=as.integer(position)-5
    query_table_clusters<-query_table[as.integer(query_table$V2) > lower_bound,]
    query_table_clusters<-query_table_clusters[as.integer(query_table_clusters$V2) < upper_bound,]
    if (length(query_table_clusters[,1]) >= 3 ){ # has to be 3 because it will always incdue itself, so if it includes itself and two others there are three SNPs within 10 kb and need to be eliminated
        SNP_clusters<-c(SNP_clusters, as.character(bowtie2_filtered_df$V1[i]))}
}
#with SNP clusters identified these should now be filtered out of the dataframe!
bowtie2_filtered_df<-bowtie2_filtered_df[!(bowtie2_filtered_df$V1 %in% SNP_clusters),]
#-------------------------------------------------------------------------------------


#------------novoalign 
contig_list<-c()
for (i in 1:length(novoalign_filtered_df[,1])){
    contig_info<-unlist(strsplit(as.character(novoalign_filtered_df$V1[i]), '_'))
    contig_identifier<-paste0(contig_info[1],'_',  contig_info[2])
    contig_list<-c(contig_list, contig_identifier)
}


#add the identifier column onto the dataframe so that we can use these names to filter out later 
novoalign_filtered_df<-cbind(novoalign_filtered_df, contig_list)




SNP_clusters<-c()
for (i in 1:length(novoalign_filtered_df[,1])){
    position<-as.integer(novoalign_filtered_df$V2[i])
    node<-novoalign_filtered_df$contig_list[i]
    query_table<-subset(novoalign_filtered_df, novoalign_filtered_df$contig_list==node)
    upper_bound=as.integer(position)+5
    lower_bound=as.integer(position)-5
    query_table_clusters<-query_table[as.integer(query_table$V2) > lower_bound,]
    query_table_clusters<-query_table_clusters[as.integer(query_table_clusters$V2) < upper_bound,]
    if (length(query_table_clusters[,1]) >= 3 ){ # has to be 3 because it will always incdue itself, so if it includes itself and two others there are three SNPs within 10 kb and need to be eliminated
        SNP_clusters<-c(SNP_clusters, as.character(novoalign_filtered_df$V1[i]))}
}
#with SNP clusters identified these should now be filtered out of the dataframe!
novoalign_filtered_df<-novoalign_filtered_df[!(novoalign_filtered_df$V1 %in% SNP_clusters),]# this line gives an error for no apparent reason 
#-------------------------------------------------------------------------------------








#------------GSNAP 
contig_list<-c()
for (i in 1:length(GSNAP_filtered_df[,1])){
    contig_info<-unlist(strsplit(as.character(GSNAP_filtered_df$V1[i]), '_'))
    contig_identifier<-paste0(contig_info[1],'_',  contig_info[2])
    contig_list<-c(contig_list, contig_identifier)
}


#add the identifier column onto the dataframe so that we can use these names to filter out later 
GSNAP_filtered_df<-cbind(GSNAP_filtered_df, contig_list)




SNP_clusters<-c()
for (i in 1:length(GSNAP_filtered_df[,1])){
    position<-as.integer(GSNAP_filtered_df$V2[i])
    node<-GSNAP_filtered_df$contig_list[i]
    query_table<-subset(GSNAP_filtered_df, GSNAP_filtered_df$contig_list==node)
    upper_bound=as.integer(position)+5
    lower_bound=as.integer(position)-5
    query_table_clusters<-query_table[as.integer(query_table$V2) > lower_bound,]
    query_table_clusters<-query_table_clusters[as.integer(query_table_clusters$V2) < upper_bound,]
    if (length(query_table_clusters[,1]) >= 3 ){ # has to be 3 because it will always incdue itself, so if it includes itself and two others there are three SNPs within 10 kb and need to be eliminated
        SNP_clusters<-c(SNP_clusters, as.character(GSNAP_filtered_df$V1[i]))}
}
#with SNP clusters identified these should now be filtered out of the dataframe!
GSNAP_filtered_df<-GSNAP_filtered_df[!(GSNAP_filtered_df$V1 %in% SNP_clusters),]
#-------------------------------------------------------------------------------------



#now add a column that describes method for each row and merge all dataframes together and then isolate duplicate values
bowtie<-c()
for (i in 1:length(bowtie2_filtered_df[,1])){
    bowtie<-c(bowtie, 'bowtie2')
}
bowtie2_filtered_df<-cbind(bowtie2_filtered_df, method=bowtie)



unique_id<-c()
for (i in 1:length(bowtie2_filtered_df[,1])){
    unique_id<-c(unique_id, paste0(bowtie2_filtered_df$contig_list[i],'_', bowtie2_filtered_df$V2[i]))
}
bowtie2_filtered_df<-cbind(bowtie2_filtered_df, idenifier=unique_id)





GSNAP<-c()
for (i in 1:length(GSNAP_filtered_df[,1])){
    GSNAP<-c(GSNAP, 'GSNAP')
}
GSNAP_filtered_df<-cbind(GSNAP_filtered_df, method=GSNAP)

unique_id<-c()
for (i in 1:length(GSNAP_filtered_df[,1])){
    unique_id<-c(unique_id, paste0(GSNAP_filtered_df$contig_list[i],'_', GSNAP_filtered_df$V2[i]))
}
GSNAP_filtered_df<-cbind(GSNAP_filtered_df, idenifier=unique_id)




novoalign<-c()
for (i in 1:length(novoalign_filtered_df[,1])){
    novoalign<-c(novoalign, 'novoalign')
}
novoalign_filtered_df<-cbind(novoalign_filtered_df, method=novoalign)


unique_id<-c()
for (i in 1:length(novoalign_filtered_df[,1])){
    unique_id<-c(unique_id, paste0(novoalign_filtered_df$contig_list[i],'_', novoalign_filtered_df$V2[i]))
}

novoalign_filtered_df<-cbind(novoalign_filtered_df, idenifier=unique_id)
#merge dataframes
merged_data_frame<-rbind(bowtie2_filtered_df,GSNAP_filtered_df,novoalign_filtered_df)

#---------------determine duplicated and triplicated SNPs 
duplicated_SNPs<-merged_data_frame$idenifier[duplicated(merged_data_frame$idenifier)]
triplicated_SNPs<-duplicated_SNPs[duplicated(duplicated_SNPs)]
duplicated_SNPs<-duplicated_SNPs[!(duplicated_SNPs %in% triplicated_SNPs)]


singlet_SNPS<-merged_data_frame$idenifier[!(merged_data_frame$idenifier %in% duplicated_SNPs | merged_data_frame$idenifier %in% triplicated_SNPs)]


merged_data_frame$V3 <-NULL
merged_data_frame$V7 <-NULL
#write out dataframes of high, medium and low quality SNPs


#---------------HQ
high_quality_data_frame<-data.frame()
for (i in triplicated_SNPs){
    high_quality_data_frame<-rbind(high_quality_data_frame, bowtie2_filtered_df[bowtie2_filtered_df$idenifier==as.character(i),])
}

for (i in triplicated_SNPs){
    high_quality_data_frame<-rbind(high_quality_data_frame, novoalign_filtered_df[novoalign_filtered_df$idenifier==i,])
}

for (i in triplicated_SNPs){
    high_quality_data_frame<-rbind(high_quality_data_frame, GSNAP_filtered_df[GSNAP_filtered_df$idenifier==i,])
}


#--------------------
#-----------------------MQ
med_quality_data_frame<-data.frame()
for (i in duplicated_SNPs){
    med_quality_data_frame<-rbind(med_quality_data_frame, bowtie2_filtered_df[bowtie2_filtered_df$idenifier==as.character(i),])
}

for (i in duplicated_SNPs){
    med_quality_data_frame<-rbind(med_quality_data_frame, novoalign_filtered_df[novoalign_filtered_df$idenifier==i,])
}

for (i in duplicated_SNPs){
    med_quality_data_frame<-rbind(med_quality_data_frame, GSNAP_filtered_df[GSNAP_filtered_df$idenifier==i,])
}

#------------------------------

#---------------LQ 
low_quality_data_frame<-data.frame()
for (i in singlet_SNPS){
    low_quality_data_frame<-rbind(low_quality_data_frame, bowtie2_filtered_df[bowtie2_filtered_df$idenifier==as.character(i),])
}

for (i in singlet_SNPS){
    low_quality_data_frame<-rbind(low_quality_data_frame, novoalign_filtered_df[novoalign_filtered_df$idenifier==i,])
}

for (i in singlet_SNPS){
    low_quality_data_frame<-rbind(low_quality_data_frame, GSNAP_filtered_df[GSNAP_filtered_df$idenifier==i,])
}


#write out SNP lists 
low_quality_data_frame$V3<-NULL
low_quality_data_frame$V7<-NULL
med_quality_data_frame$V3<-NULL
med_quality_data_frame$V7<-NULL
high_quality_data_frame$V3<-NULL
high_quality_data_frame$V7<-NULL

col_names<-c('contig','position','ref','alternate','qual','info','format','ploidy_info','contig', 'mapper','unique_SNP_ID')
colnames(low_quality_data_frame)<-col_names
colnames(med_quality_data_frame)<-col_names
colnames(high_quality_data_frame)<-col_names
write.csv(high_quality_data_frame, 'high_quality_SNPS.csv')
write.csv(med_quality_data_frame, 'med_quality_SNPS.csv')
write.csv(low_quality_data_frame, 'low_quality_SNPS.csv')






