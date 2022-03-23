#March 24, 2022
#9k snps common 1000G and Jomon _Afanasievo
#genotypes extracted
#polygenic score calculated
#EA proxy replication for all snps...
#uncomment to save plots etc
#Author kaisar dauyey

rm(list = ls())
library(dplyr)
library("ggsci")
library(gridExtra)
library(ggplot2)
library(cowplot)

library(tidyverse)
library(ggpubr)
library(rstatix)

#setwd("data/")
setwd("/Users/kaisar_dauyey/testdir/savage/plink_intelligence/scripts/paper/data")

#load 1000gp phase 3 list of people
phase3<- read.table("phase3_people.csv", 
                    sep = "\t" , header = F,
                    na.strings ="", stringsAsFactors= F, 
)

#merged PCA
#1000G_jomon_12k_AFA
#1402 SNPS and 9k SNPS



# read in data 1402
pca_1402<- read_table2("plink/merged_update_1000G_J_Afa_1880.eigenvec", col_names = FALSE)
eigenval_1402 <- scan("plink/merged_update_1000G_J_Afa_1880.eigenval")

# read in data 9k
pca_9k<- read_table2("plink/merged_update_1000G_J_Afa.eigenvec", col_names = FALSE)
eigenval_9k <- scan("plink/merged_update_1000G_J_Afa.eigenval")



# sort out the pca data
# remove nuisance column
pca_1402 <- pca_1402[,-1]
names(pca_1402)[1] <- "ind"
# set names per phase3


names(pca_1402)[2:ncol(pca_1402)] <- paste0("PC", 1:(ncol(pca_1402)-1))
pops<-data.frame(phase3$V2,phase3$V3,stringsAsFactors=FALSE)
extra1<-c('Jomon', "Jomon")
extra2<-c('I3388.Mother', "Afanasievo_m")
extra3<-c('I3949.Son1', "Afanasievo_s1")
extra4<-c('I3950.Father', "Afanasievo_f")
extra5<-c('I6714.Son2', "Afanasievo_s2")



pops<-rbind(pops,extra1,extra2,extra3,extra4,extra5)
#pops<-rbind(pops,extra1)

tail(pops)
names(pops)<-c("ind","popa")
pops<-data.frame(pops)
m1_1402 <- merge(pca_1402, pops, by.x = "ind",all.x = TRUE)
pca_1402<-m1_1402

tail(pca_1402)

#find Jomon
which(pca_1402$popa == "Jomon")
#1661
which(pca_1402$popa == "Afanasievo_s1")
#1658
which(pca_1402$popa == "Afanasievo_s2")
#1660
which(pca_1402$popa == "Afanasievo_m")
#1657
which(pca_1402$popa == "Afanasievo_f")
#1659

# sort out the individual species and pops
# spp
spp_1402 <- rep(NA, length(pca_1402$popa))

#pca[,1]<- phase3[,3]
spp_1402[grep("CHB", pca_1402$popa)] <- "EAS"
spp_1402[grep("JPT", pca_1402$popa)] <- "EAS"
spp_1402[grep("CHS", pca_1402$popa)] <- "EAS"
spp_1402[grep("CDX", pca_1402$popa)] <- "EAS"
spp_1402[grep("KHV", pca_1402$popa)] <- "EAS"
spp_1402[grep("CEU", pca_1402$popa)] <- "EUR"
spp_1402[grep("TSI", pca_1402$popa)] <- "EUR"
spp_1402[grep("FIN", pca_1402$popa)] <- "EUR"
spp_1402[grep("GBR", pca_1402$popa)] <- "EUR"
spp_1402[grep("IBS", pca_1402$popa)] <- "EUR"
spp_1402[grep("YRI", pca_1402$popa)] <- "AFR"
spp_1402[grep("LWK", pca_1402$popa)] <- "AFR"
spp_1402[grep("GWD", pca_1402$popa)] <- "AFR"
spp_1402[grep("MSL", pca_1402$popa)] <- "AFR"
spp_1402[grep("ESN", pca_1402$popa)] <- "AFR"
spp_1402[grep("ASW", pca_1402$popa)] <- "AFR"
spp_1402[grep("ACB", pca_1402$popa)] <- "AFR"
spp_1402[grep("MXL", pca_1402$popa)] <- "AMR"
spp_1402[grep("PUR", pca_1402$popa)] <- "AMR"
spp_1402[grep("CLM", pca_1402$popa)] <- "AMR"
spp_1402[grep("PEL", pca_1402$popa)] <- "AMR"
spp_1402[grep("GIH", pca_1402$popa)] <- "SAS"
spp_1402[grep("PJL", pca_1402$popa)] <- "SAS"
spp_1402[grep("BEB", pca_1402$popa)] <- "SAS"
spp_1402[grep("STU", pca_1402$popa)] <- "SAS"
spp_1402[grep("ITU", pca_1402$popa)] <- "SAS"
spp_1402[grep("Jomon", pca_1402$popa)] <- "F23"
spp_1402[grep("Afanasievo_s2", pca_1402$popa)] <- "AS2"
spp_1402[grep("Afanasievo_s1", pca_1402$popa)] <- "AS1"
spp_1402[grep("Afanasievo_m", pca_1402$popa)] <- "AM"
spp_1402[grep("Afanasievo_f", pca_1402$popa)] <- "AF"

# remake data.frame
pca_1402 <- as_tibble(data.frame(pca_1402, spp_1402))
tail(pca_1402)
pca_1402$spp_1402 <- factor(pca_1402$spp_1402)

# first convert to percentage variance explained
pve_1402 <- data.frame(PC = 1:20, pve = eigenval_1402/sum(eigenval_1402)*100)
# make plot
a_1402 <- ggplot(pve_1402, aes(PC, pve)) + geom_bar(stat = "identity")
a_1402 + ylab("Percentage variance explained") + theme_light()
ala_1402<-a_1402 + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve_1402$pve)
#ggsave("pve4.png", ala_1402, width=7, height=5, units="in", dpi=500)

#new pca



#1000G_jomon_ 9k_AFA
# sort out the pca data
# remove nuisance column
pca_9k <- pca_9k[,-1]
names(pca_9k)[1] <- "ind"
# set names per phase3


names(pca_9k)[2:ncol(pca_9k)] <- paste0("PC", 1:(ncol(pca_9k)-1))
pops<-data.frame(phase3$V2,phase3$V3,stringsAsFactors=FALSE)
extra1<-c('Jomon', "Jomon")
extra2<-c('I3388.Mother', "Afanasievo_m")
extra3<-c('I3949.Son1', "Afanasievo_s1")
extra4<-c('I3950.Father', "Afanasievo_f")
extra5<-c('I6714.Son2', "Afanasievo_s2")



pops<-rbind(pops,extra1,extra2,extra3,extra4,extra5)
#pops<-rbind(pops,extra1)

tail(pops)
names(pops)<-c("ind","popa")
pops<-data.frame(pops)
m1_9k <- merge(pca_9k, pops, by.x = "ind",all.x = TRUE)
pca_9k<-m1_9k

tail(pca_9k)

#find Jomon
which(pca_9k$popa == "Jomon")
#1661
which(pca_9k$popa == "Afanasievo_s1")
#1658
which(pca_9k$popa == "Afanasievo_s2")
#1660
which(pca_9k$popa == "Afanasievo_m")
#1657
which(pca_9k$popa == "Afanasievo_f")
#1659

# sort out the individual species and pops
# spp
spp_9k <- rep(NA, length(pca_9k$popa))

#pca[,1]<- phase3[,3]
spp_9k[grep("CHB", pca_9k$popa)] <- "EAS"
spp_9k[grep("JPT", pca_9k$popa)] <- "EAS"
spp_9k[grep("CHS", pca_9k$popa)] <- "EAS"
spp_9k[grep("CDX", pca_9k$popa)] <- "EAS"
spp_9k[grep("KHV", pca_9k$popa)] <- "EAS"
spp_9k[grep("CEU", pca_9k$popa)] <- "EUR"
spp_9k[grep("TSI", pca_9k$popa)] <- "EUR"
spp_9k[grep("FIN", pca_9k$popa)] <- "EUR"
spp_9k[grep("GBR", pca_9k$popa)] <- "EUR"
spp_9k[grep("IBS", pca_9k$popa)] <- "EUR"
spp_9k[grep("YRI", pca_9k$popa)] <- "AFR"
spp_9k[grep("LWK", pca_9k$popa)] <- "AFR"
spp_9k[grep("GWD", pca_9k$popa)] <- "AFR"
spp_9k[grep("MSL", pca_9k$popa)] <- "AFR"
spp_9k[grep("ESN", pca_9k$popa)] <- "AFR"
spp_9k[grep("ASW", pca_9k$popa)] <- "AFR"
spp_9k[grep("ACB", pca_9k$popa)] <- "AFR"
spp_9k[grep("MXL", pca_9k$popa)] <- "AMR"
spp_9k[grep("PUR", pca_9k$popa)] <- "AMR"
spp_9k[grep("CLM", pca_9k$popa)] <- "AMR"
spp_9k[grep("PEL", pca_9k$popa)] <- "AMR"
spp_9k[grep("GIH", pca_9k$popa)] <- "SAS"
spp_9k[grep("PJL", pca_9k$popa)] <- "SAS"
spp_9k[grep("BEB", pca_9k$popa)] <- "SAS"
spp_9k[grep("STU", pca_9k$popa)] <- "SAS"
spp_9k[grep("ITU", pca_9k$popa)] <- "SAS"
spp_9k[grep("Jomon", pca_9k$popa)] <- "F23"
spp_9k[grep("Afanasievo_s2", pca_9k$popa)] <- "AS2"
spp_9k[grep("Afanasievo_s1", pca_9k$popa)] <- "AS1"
spp_9k[grep("Afanasievo_m", pca_9k$popa)] <- "AM"
spp_9k[grep("Afanasievo_f", pca_9k$popa)] <- "AF"

# remake data.frame
pca_9k <- as_tibble(data.frame(pca_9k, spp_9k))
tail(pca_9k)
pca_9k$spp_9k <- factor(pca_9k$spp_9k)

# first convert to percentage variance explained
pve_9k <- data.frame(PC = 1:20, pve = eigenval_9k/sum(eigenval_9k)*100)
# make plot
a_9k <- ggplot(pve_9k, aes(PC, pve)) + geom_bar(stat = "identity")
a_9k + ylab("Percentage variance explained") + theme_light()
ala_9k<-a_9k + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve_9k$pve)
#ggsave("pve4.png", ala_9k, width=7, height=5, units="in", dpi=500)



#plotting
#ka_1402<- c(15+ nlevels(pca_1402$spp_1402))

# plot pca
b_1402 <- ggplot(pca_1402, aes(PC1, PC2, col = spp_1402 )) +
  geom_point(size = 0.5) +
  geom_point(data=pca_1402[1661, ], size=3,) +
  geom_point(data=pca_1402[1658, ], size=3,) +
  geom_point(data=pca_1402[1657, ], size=3,) +
  geom_point(data=pca_1402[1660, ], size=3,) +
  geom_point(data=pca_1402[1659, ], size=3,) +
  scale_color_npg()  

b_1402 <- b_1402 + coord_equal()+ theme_cowplot(font_size = 7, line_size = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
  


b_1402 
aca_1402<-b_1402 

#9k trial
#ka_9k<- c(15+ nlevels(pca_9k$spp_9k))

# plot pca
# plot pca
b_9k <- ggplot(pca_9k, aes(PC1, PC2, col = spp_9k )) +
  geom_point(size = 0.5) +
  geom_point(data=pca_9k[1661, ], size=3) +
  geom_point(data=pca_9k[1658, ], size=3) +
  geom_point(data=pca_9k[1657, ], size=3) +
  geom_point(data=pca_9k[1660, ], size=3) +
  geom_point(data=pca_9k[1659, ], size=3) +
  scale_color_npg()  

b_9k <- b_9k + coord_equal()+ theme_cowplot(font_size = 7, line_size = 1)+
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


b_9k 
aca_9k<-b_9k

Figure2 <- plot_grid(aca_1402, aca_9k, labels = "AUTO", ncol = 1, align = 'v')

Figure2

#ggsave("~/Desktop/Figure2.pdf", Figure2, width=8.5, height=13, units="cm", dpi=300)
