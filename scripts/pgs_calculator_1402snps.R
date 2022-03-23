#March 24, 2022
#intelligence_counter_with_Beta_score_
#1000G with Jomon Afanasievo
#1402 snps used total
#positive and negative snp count included
#Author Kaisar Dauyey



rm(list = ls())
library(dplyr)
library(plyr)
library(ggplot2)

setwd("/Ancient_intelligence/")

#load 1000gp phase 3 list of people
phase3<- read.table("data/phase3_people.csv", 
                    sep = "\t" , header = F,
                    na.strings ="", stringsAsFactors= F, 
)
head(phase3$V3)
#REFERENCE load from SAVAGE et.al

ref12k<- read.table("data/total_12k_snp_alleles_beta_uk_EA.csv", 
                    sep = "," , header = T,
                    na.strings ="", stringsAsFactors= F, 
)
ref<-data.frame(ref12k)
head(ref)



ref$Effect_allele = toupper(ref$Effect_allele)

# load cleaned up ped and map files with 12k  - unmapped 26 snps abscent in ancient individuals


G_J_12k.ped<- read.csv("data/plink/merged_update_1000G_J_Afa_1880.ped", sep = ' ', header=FALSE, stringsAsFactors=FALSE, 
                       colClasses = c("character"))

mapa_ALL_1000GP <- read.csv("data/plink/merged_update_1000G_J_Afa_1880.map", sep = '', header=FALSE, stringsAsFactors=FALSE, 
                            colClasses = c("character"))




G_J_12k.ped<-G_J_12k.ped[,-c(2:6)]

str(G_J_12k.ped)

#check for NAs


# order 1000GP check - take populations

#check missing SNPS in 12k out of 12k (EA total ---- 1878 subset by plink)
index1<- which(!ref$BP %in% mapa_ALL_1000GP[,4] )
#index2<- which(!mapa_ALL_1000GP[,4] %in% ref$BP)

#gomi1 <-as.data.frame(ref[3][index1,])
#gomi2 <-as.data.frame(mapa_ALL_1000GP[2][index2,])


ref2<- ref[-index1,]
#mapa_ALL_1000GP2<-mapa_ALL_1000GP[-index2,]

#clean_mapa<-mapa_ALL_1000GP2[2]
#write.csv(clean_mapa, "cleaning_gomi.csv")

head(ref2$SNP)

colnames(mapa_ALL_1000GP)[2] <- "SNP"

table_six <- inner_join(x = mapa_ALL_1000GP, 
                        y = ref2,
                        by = "SNP") %>%
  unique()


str(table_six)

names<-c("subject",rep(table_six$SNP, each=2))
head(names)

colnames(G_J_12k.ped)<-names


#counting.with effect...

#counting.with effect...
m1<-data.frame(G_J_12k.ped)
#head(m1)
pops<-data.frame(phase3$V2,phase3$V3,stringsAsFactors=FALSE)
extra1<-c('Jomon', "Jomon")
extra2<-c('I3388.Mother', "Afanasievo_m")
extra3<-c('I3949.Son1', "Afanasievo_s1")
extra4<-c('I3950.Father', "Afanasievo_f")
extra5<-c('I6714.Son2', "Afanasievo_s2")

pops<-rbind(pops,extra1,extra2,extra3,extra4,extra5)

names(pops)<-c("subject","popa")
m2 <- merge(m1, pops, by.x = "subject",all.x = TRUE)

df <- data.frame(matrix(ncol = 7, nrow = nrow(m2)))
x <- c("subject", "total_ALL_1000GP","positive_count","negative_count", "positive_effect",
       "negative_effect", "effect_sum")
colnames(df) <- x

#remove last column of m2 - population names for calculation
df[,1]<- m2[,2806]

analyze3 <- function(file1) {
  
# ALSO remove last column of m2 - population names for calculation
#change end product CSV name of calculated results   
  df2 <- m2[, -c(1,2806)]
  ref_test1<- as.matrix(as.character(rep(table_six$Effect_allele,each=2)))
  ref_effect1<- as.matrix(as.numeric(rep(table_six$Beta_EA,each=2)))
  #ref_effect_freq<- as.matrix(as.numeric(rep(ref2$EAF_HRC,each=2)))
  #no EAF_HRC used for score calculation
  for (v in 1:nrow(df2)){
    
    trial1<-as.matrix(as.character(df2[v,]))
    
    picker1<-which(trial1 == ref_test1)
    effect_subject1<-ref_effect1[picker1]
    
    #effect_subject1<-ref_effect1[picker1]*ref_effect_freq[picker1]
    
    
    snp_count1<- length(picker1)
    
    effector1<-sum(effect_subject1)
    
    
    positive_count<-length(effect_subject1[which(effect_subject1>0)])
    negative_count<-length(effect_subject1[which(effect_subject1<0)])
    positive_effect<-sum(effect_subject1[which(effect_subject1>0)])
    negative_effect<-sum(effect_subject1[which(effect_subject1<0)])
    
    df[v,2] <- as.numeric(snp_count1)
    df[v,3] <- as.numeric(positive_count)
    df[v,4] <- as.numeric(negative_count)
    df[v,5] <- as.numeric(positive_effect)
    df[v,6] <- as.numeric(negative_effect)
    df[v,7] <- effector1
  }
  write.csv(  df,'PGS_1000G_J_Afa_1402.csv',row.names = FALSE)
  book1<<-data.frame(df)
  
}

analyze3(m1)

hist(book1$effect_sum)

#exploring the ancient genomes PGS values for 1402 snps

which(book1$subject == "Jomon")

book1[1661,]

book1[1657,]
book1[1658,]
book1[1659,]
book1[1660,]

