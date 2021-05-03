#1check ancestral state of snps reported by Savage by ORTHEUS package
#May 3 , 2021
#11 627 out 12 037 common in Jomon and 1000 Genome subjects were tested
# skip to "plot with data calculated" to avoid derivation of ancestry

#Author kaisar dauyey

library(biomaRt)
library(ggplot2)
library(dplyr)
library(viridis)
setwd("/Users/kaisar_dauyey/testdir/savage/")
library(janitor)


#Load reference 

human_variation = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp") 
#GRCh37/hg19 as in Savage et al 2018
#grch37 = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
#read data table with 12 037 snps in total
data <- read.csv("12k_snp_alleles.csv")
head(data)
dim(data)
snp_ids<-as.vector(data$SNP)

#put attributes
snp_attributes = c("refsnp_id", 
                   "chr_name", 
                   "chrom_start", 
                   "chrom_end",
                   "allele",
                   "allele_1");

#use asian server
asia_mart <- useEnsembl(biomart = "ensembl", mirror = "asia")
snp_locations = getBM(attributes = snp_attributes, 
                      filters = "snp_filter", 
                      values = snp_ids, 
                      mart = human_variation,
                      );

#remove duplicates - remove them as they contain characters in chr name

value <-  as.numeric( as.character(snp_locations$chr_name) )
noninteger <- value %% 1 != 0   # see if there's a fractional part
noninteger <- noninteger | is.na(noninteger)  # get rid of the NA's
data3<- snp_locations[!noninteger,]

#order data based on snp names
data4<-data3[order(data3$refsnp_id),]
head(data4)

#vector_anc <- as.character(data_anc$allele_1)
#structure(vector_anc)
empty_anc<-which(data4$allele_1 == '')  
#extract to csv

write.csv(data4, "snps_11627_ancestry_hg38.csv")

#plot with data calculated

data_ref<- read.csv("12k_snp_alleles.csv")
data_anc<- read.csv("snps_11627_ancestry_hg38.csv")
head(data_ref)
head(data_anc)


data_ref$Effect_allele = toupper(data_ref$Effect_allele)

#no need to allign - do merge by matching SNPs
is.null(data_anc)

#merge by match...
data_merged<-merge(x = data_ref, y = data_anc,
                   by.x = "SNP", by.y = "refsnp_id")
head(data_merged)
dim(data_merged)
#checkNA
sum(is.na(data_merged$allele_1))

#all alleles
sum_all<-data.frame(table(data_merged$Effect_allele == data_merged$allele_1))
sum_all[1] <- lapply(sum_all[1], as.character)
sum_all <- rbind(sum_all, c(sum(sum_all$Freq)))


d1<-c("Derived", "Ancestral", "Total")
sum_all<-cbind(sum_all,d1)
head(sum_all)


#subset all +
data_pos<-data_merged[data_merged$Zscore>0,]
head(data_pos)
#subset all -
data_neg<-data_merged[!data_merged$Zscore>0,]
head(data_neg)


#counter +
sum_pos<-as.data.frame(table(data_pos$Effect_allele == data_pos$allele_1))
sum_pos[1] <- lapply(sum_pos[1], as.character)
sum_pos <- rbind(sum_pos, c(sum(sum_pos$Freq)))

sum_pos<-cbind(sum_pos,d1)
head(sum_pos)

#counter -
sum_neg<-as.data.frame(table(data_neg$Effect_allele == data_neg$allele_1))
sum_neg[1] <- lapply(sum_neg[1], as.character)
sum_neg <- rbind(sum_neg, c(sum(sum_neg$Freq)))

sum_neg<-cbind(sum_neg,d1)
head(sum_neg)



counts<-bind_rows(sum_all,sum_pos,sum_neg)
dim_desc(counts)
counts



allele <- c(rep("Absolute count" , 3) , rep("Positive" , 3) , rep("Negative" , 3) )
data_test <- data.frame(counts,allele)
head(data_test)

#write ancestry statistics
#write.csv(data_test,"ancestry_12k.csv")



