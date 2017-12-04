#====================Take user input VCFs
args = commandArgs(trailingOnly=TRUE)
dir=args[1]
reference=args[2]
#dir='/Users/Harley/Desktop/SNP_analysis/gzhou_3'
#reference='/Users/Harley/Desktop/SNP_analysis/gzhou_2/LP_Philly.fasta'
library(seqinr)
A<-read.fasta(reference, as.string = T, seqonly = T, strip.desc = T)
length<-as.integer(length(unlist(strsplit(as.character(A[[1]][1]),''))))
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
#check if the dataframe is empty
if (!(is.null(count.fields(files[1])))){ # this line checks to see if the vcf is empty, if it isnt we populate the dataframe
    bowtie2_filtered_df<-data.frame()
    bowtie2_df<-read.table(files[1]) 
    for (i in 1:length(bowtie2_df[,1])){
        coverage<-as.integer(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[6])
        identifier<-cat(as.character(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[1]),'_',as.character(unlist(strsplit(as.character(bowtie2_df$V1[i]),'_'))[2]),'_',as.character(bowtie2_df$V2[i]))
        position<-as.integer(bowtie2_df$V2[i])
        if(position  > distance_to_end_of_contig & (length-position > distance_to_end_of_contig)){
            bowtie2_filtered_df<-rbind(bowtie2_filtered_df, bowtie2_df[i,])
        }
        
    } 
}





#--------------------------------------------------------------
#------------------------------GSNAP
if (!(is.null(count.fields(files[2])))){
    GSNAP_df<-read.table(files[2])
    GSNAP_filtered_df<-data.frame()
    for (i in 1:length(GSNAP_df[,1])){
        length<-23091291
        coverage<-as.integer(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[6])
        identifier<-cat(as.character(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[1]),'_',as.character(unlist(strsplit(as.character(GSNAP_df$V1[i]),'_'))[2]),'_',as.character(GSNAP_df$V2[i]))
        position<-as.integer(GSNAP_df$V2[i])
        if(position  > distance_to_end_of_contig & (length-position > distance_to_end_of_contig)){
            GSNAP_filtered_df<-rbind(GSNAP_filtered_df, GSNAP_df[i,])
        }
        
    }
}
#-----------------------------------------------------------------
#--------------------------------novoalign
if (!(is.null(count.fields(files[3])))){
    novoalign_df<-read.table(files[3])
    novoalign_filtered_df<-data.frame()
    for (i in 1:length(novoalign_df[,1])){
        length<-2301291
        coverage<-as.integer(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[6])
        identifier<-cat(as.character(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[1]),'_',as.character(unlist(strsplit(as.character(novoalign_df$V1[i]),'_'))[2]),'_',as.character(novoalign_df$V2[i]))
        position<-as.integer(novoalign_df$V2[i])
        if(position  > distance_to_end_of_contig & (length-position > distance_to_end_of_contig)){
            novoalign_filtered_df<-rbind(novoalign_filtered_df, novoalign_df[i,])
        }
        
    }}
#-------------------------------------------------------------------

#-------------Now want to filter SNPs found in a cluster, this will be loosely defined as more than 3 SNPs within a 10 bp window 
#the filtration for clusters comes from here : Ren, C. et al. SNP discovery by high-throughput sequencing in soybean. BMC Genomics 11, 469 (2010).

#get the identifier --> node and position for all SNPs and then filter within tables first to remove stretches of SNPs 
#then split SNPs into high quality (found in all three) or  medium qualtiy (found in 2/3) or low quality (found in 1/3) samples
#---------------------------------------------SNP cluster program for finding clusters of SNPs and excluding them (likely mapping errors)

#------------Bowtie 
#get bowtie SNPs
if (exists('bowtie2_filtered_df')){ #added this if exists check to see if the dataframe is actually present before I start doing stuff with it 
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
    
}
#-------------------------------------------------------------------------------------
#------------novoalign 
if (exists('novoalign_filtered_df')){
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
}
#-------------------------------------------------------------------------------------
#------------GSNAP 
if (exists('GSNAP_filtered_df')){
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
}
#-------------------------------------------------------------------------------------
#now add a column that describes method for each row and merge all dataframes together and then isolate duplicate values
#note that in V1 this was its own section, but I am moving it into the exists statements so it doesnt make changes unless the dataframes exist 
#bowtie<-c()
#for (i in 1:length(bowtie2_filtered_df[,1])){
#    bowtie<-c(bowtie, 'bowtie2')
#}
#bowtie2_filtered_df<-cbind(bowtie2_filtered_df, method=bowtie)

#unique_id<-c()
#for (i in 1:length(bowtie2_filtered_df[,1])){
#    unique_id<-c(unique_id, paste0(bowtie2_filtered_df$contig_list[i],'_', bowtie2_filtered_df$V2[i]))
#}
#bowtie2_filtered_df<-cbind(bowtie2_filtered_df, idenifier=unique_id)


#----------------------------------------------
#GSNAP<-c()
#for (i in 1:length(GSNAP_filtered_df[,1])){
#    GSNAP<-c(GSNAP, 'GSNAP')
#}
#GSNAP_filtered_df<-cbind(GSNAP_filtered_df, method=GSNAP)
#
#unique_id<-c()
#for (i in 1:length(GSNAP_filtered_df[,1])){
#    unique_id<-c(unique_id, paste0(GSNAP_filtered_df$contig_list[i],'_', GSNAP_filtered_df$V2[i]))
#}
#GSNAP_filtered_df<-cbind(GSNAP_filtered_df, idenifier=unique_id)
#---------------------------------------- 

#-----------------------------------------
#novoalign<-c()
#for (i in 1:length(novoalign_filtered_df[,1])){
#    novoalign<-c(novoalign, 'novoalign')
#}
#novoalign_filtered_df<-cbind(novoalign_filtered_df, method=novoalign)
##
#
#unique_id<-c()
#for (i in 1:length(novoalign_filtered_df[,1])){
#    unique_id<-c(unique_id, paste0(novoalign_filtered_df$contig_list[i],'_', novoalign_filtered_df$V2[i]))
#}
#
#novoalign_filtered_df<-cbind(novoalign_filtered_df, idenifier=unique_id)
#-----------------------------------------
#merge dataframes


#---------------UPDATING NOTE: Now dataframes that have data are created. Need to merge dataframes but only if they exist
#prior to merging the dataframes into a dataframe need to check if they exist
#make a operator like bowtie_exists if true its merged to DF if false its not

bowtie_exists<-exists('bowtie2_filtered_df')
GSNAP_exists<-exists('GSNAP_filtered_df')
novoalign_exists<-exists('novoalign_filtered_df')
#------ now that we know which datasets exist we can merge them together if they exist
#can do this stepwise to avoid confusion
merged_data_frame<-data.frame()
if (bowtie_exists){
    merged_data_frame<-rbind(merged_data_frame, bowtie2_filtered_df)
}

if (GSNAP_exists){
    merged_data_frame<-rbind(merged_data_frame, GSNAP_filtered_df)
}

if (novoalign_exists){
    merged_data_frame<-rbind(merged_data_frame, novoalign_filtered_df)
}



#---------------determine duplicated and triplicated SNPs 
duplicated_SNPs<-merged_data_frame$idenifier[duplicated(merged_data_frame$idenifier)]
triplicated_SNPs<-duplicated_SNPs[duplicated(duplicated_SNPs)]#fun trick, the doubly duplciated's are triplicates thus in all three mappers 
duplicated_SNPs<-duplicated_SNPs[!(duplicated_SNPs %in% triplicated_SNPs)]#those not triplicated are just duplicated


singlet_SNPS<-merged_data_frame$idenifier[!(merged_data_frame$idenifier %in% duplicated_SNPs | merged_data_frame$idenifier %in% triplicated_SNPs)]#those not dup or trip are singlets and Low quality


merged_data_frame$V3 <-NULL# remove . columns
merged_data_frame$V7 <-NULL # remove . columns
#write out dataframes of high, medium and low quality SNPs





#UPDATED NOTES: BECAUSE I CHANGED THE CODE TO CHECK IF A DATAFRAME EXISTS I ALSO NEED TO APPEND THINGS 
#INDIVIDUALLY TO THE HIGH QUALITY DATA FRAME LIKE ABOVE IF EACH THING EXISTS, IF NOT DONT APPEND
#---------------HQ
# WE ALREADY KNOW IF THE DATAFRAME EXISTS SO WE CAN JUST MODIFY THIS CODE 
high_quality_data_frame<-data.frame()
if (bowtie_exists){
    for (i in triplicated_SNPs){
        high_quality_data_frame<-rbind(high_quality_data_frame, bowtie2_filtered_df[bowtie2_filtered_df$idenifier==as.character(i),])
        
    }
}

if (novoalign_exists){
    for (i in triplicated_SNPs){
        high_quality_data_frame<-rbind(high_quality_data_frame, novoalign_filtered_df[novoalign_filtered_df$idenifier==i,])
    }
}
if (GSNAP_exists){
    for (i in triplicated_SNPs){
        high_quality_data_frame<-rbind(high_quality_data_frame, GSNAP_filtered_df[GSNAP_filtered_df$idenifier==i,])
    }
}


#--------------------
#-----------------------MQ
med_quality_data_frame<-data.frame()
if (bowtie_exists){
    for (i in duplicated_SNPs){
        med_quality_data_frame<-rbind(med_quality_data_frame, bowtie2_filtered_df[bowtie2_filtered_df$idenifier==as.character(i),])
    }
}

if (novoalign_exists){
    for (i in duplicated_SNPs){
        med_quality_data_frame<-rbind(med_quality_data_frame, novoalign_filtered_df[novoalign_filtered_df$idenifier==i,])
    }
}
if (GSNAP_exists){
    for (i in duplicated_SNPs){
        med_quality_data_frame<-rbind(med_quality_data_frame, GSNAP_filtered_df[GSNAP_filtered_df$idenifier==i,])
    }
}
#------------------------------
#---------------LQ 
low_quality_data_frame<-data.frame()
if (bowtie_exists){
    for (i in singlet_SNPS){
        low_quality_data_frame<-rbind(low_quality_data_frame, bowtie2_filtered_df[bowtie2_filtered_df$idenifier==as.character(i),])
    }
}
if (novoalign_exists){
    for (i in singlet_SNPS){
        low_quality_data_frame<-rbind(low_quality_data_frame, novoalign_filtered_df[novoalign_filtered_df$idenifier==i,])
    }
}
if (GSNAP_exists){
    for (i in singlet_SNPS){
        low_quality_data_frame<-rbind(low_quality_data_frame, GSNAP_filtered_df[GSNAP_filtered_df$idenifier==i,])
    }
}

#write out SNP lists 
#now I should check whether a output file exists, if it does delete the columns, add a header, and write to file
sample_name_global<-unlist(strsplit(sample_names[1],'-'))[1]


#if the low quality dataframe exists, write it to file 
if (exists('low_quality_data_frame')& length(low_quality_data_frame) != 0){
    low_quality_data_frame$V3<-NULL
    low_quality_data_frame$V7<-NULL
    col_names<-c('contig','position','ref','alternate','qual','info','format','ploidy_info','contig', 'mapper','unique_SNP_ID')
    colnames(low_quality_data_frame)<-col_names 
    write.csv(low_quality_data_frame, paste0('low_quality_SNPS','_',sample_name_global,'.csv'))
}

#if the med quality dataframe exists, write it to file 
if (exists('med_quality_data_frame') & length(med_quality_data_frame) != 0){
    med_quality_data_frame$V3<-NULL
    med_quality_data_frame$V7<-NULL
    col_names<-c('contig','position','ref','alternate','qual','info','format','ploidy_info','contig', 'mapper','unique_SNP_ID')
    colnames(med_quality_data_frame)<-col_names
    write.csv(med_quality_data_frame, paste0('med_quality_SNPS','_',sample_name_global,'.csv'))
}

#if the high quality dataframe exists, write it to file 
if (exists('high_quality_data_frame') & length(high_quality_data_frame) != 0){
    high_quality_data_frame$V3<-NULL
    high_quality_data_frame$V7<-NULL
    col_names<-c('contig','position','ref','alternate','qual','info','format','ploidy_info','contig', 'mapper','unique_SNP_ID')
    colnames(high_quality_data_frame)<-col_names
    write.csv(high_quality_data_frame, paste0('high_quality_SNPS','_',sample_name_global,'.csv'))
}





