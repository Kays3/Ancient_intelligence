#PCA_nean_deni_jomon_1000GP_intell
#Chimpanzee included as reference
#neanderthal and denisovan variants alligne SNPS 72 used
#October 22, 2020
#Author kaisar dauyey



rm(list = ls())
library(dplyr)
library(ggplot2)
library(adegenet)
library(tidyverse)
library(plyr)
library(readr)




setwd("~/pca")


#load 1000gp phase 3 list of people
phase3<- read.table("phase3_people.csv", 
                    sep = "\t" , header = F,
                    na.strings ="", stringsAsFactors= F, 
)

#merged PCA from PLINK analysis is available 1000G_jomon_nean_chimp


# read in data from Altai Neanderthal, Denisovan, Jomon, 1000 Genomes 
pca <- read_table2("./72_snps/1000G_jomon_nean_chimp.eigenvec", col_names = FALSE)
eigenval <- scan("./72_snps/1000G_jomon_nean_chimp.eigenval")


# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
names(pca)[1] <- "ind"
# set names per phase3


names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pops<-data.frame(phase3$V2,phase3$V3,stringsAsFactors=FALSE)
extra1<-c('AltaiNea', "AltaiNea")
extra2<-c('Denisovan', "Denisovan")
extra3<-c('Jomon', "Jomon")
extra4<-c('taweh', "taweh")



pops<-rbind(pops,extra1,extra2,extra3,extra4)
tail(pops)
names(pops)<-c("ind","popa")
pops<-data.frame(pops)
m1 <- merge(pca, pops, by.x = "ind",all.x = TRUE)
pca<-m1

tail(pca)
which(pca$popa == "Jomon")
#1659
which(pca$popa == "AltaiNea")
#1
which(pca$popa == "Denisovan")
#2
which(pca$popa == "taweh")
#2508


# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$popa))

#pca[,1]<- phase3[,3]
spp[grep("CHB", pca$popa)] <- "East_Asia"
spp[grep("JPT", pca$popa)] <- "East_Asia"
spp[grep("CHS", pca$popa)] <- "East_Asia"
spp[grep("CDX", pca$popa)] <- "East_Asia"
spp[grep("KHV", pca$popa)] <- "East_Asia"
spp[grep("CEU", pca$popa)] <- "Europe"
spp[grep("TSI", pca$popa)] <- "Europe"
spp[grep("FIN", pca$popa)] <- "Europe"
spp[grep("GBR", pca$popa)] <- "Europe"
spp[grep("IBS", pca$popa)] <- "Europe"
spp[grep("YRI", pca$popa)] <- "Africa"
spp[grep("LWK", pca$popa)] <- "Africa"
spp[grep("GWD", pca$popa)] <- "Africa"
spp[grep("MSL", pca$popa)] <- "Africa"
spp[grep("ESN", pca$popa)] <- "Africa"
spp[grep("ASW", pca$popa)] <- "Africa"
spp[grep("ACB", pca$popa)] <- "Africa"
spp[grep("MXL", pca$popa)] <- "Ad_Mixed_American"
spp[grep("PUR", pca$popa)] <- "Ad_Mixed_American"
spp[grep("CLM", pca$popa)] <- "Ad_Mixed_American"
spp[grep("PEL", pca$popa)] <- "Ad_Mixed_American"
spp[grep("GIH", pca$popa)] <- "South_Asia"
spp[grep("PJL", pca$popa)] <- "South_Asia"
spp[grep("BEB", pca$popa)] <- "South_Asia"
spp[grep("STU", pca$popa)] <- "South_Asia"
spp[grep("ITU", pca$popa)] <- "South_Asia"
spp[grep("AltaiNea", pca$popa)] <- "AltaiNea"
spp[grep("Denisovan", pca$popa)] <- "Denisovan"
spp[grep("Jomon", pca$popa)] <- "Jomon"
spp[grep("taweh", pca$popa)] <- "Chimpanzee"



# remake data.frame
pca <- as_tibble(data.frame(pca, spp))
tail(pca)
pca$spp <- factor(pca$spp)


# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
ala<-a + ylab("Percentage variance explained") + theme_light()
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
#ggsave("pve4.png", ala, width=7, height=5, units="in", dpi=500)

ka<- c(13+ nlevels(pca$spp))

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = spp )) +
  geom_point(size = 2,alpha = 1/10) +
  geom_point(data=pca[1659, ], size=5,alpha = 9/10) +
  geom_point(data=pca[2508, ], size=5,alpha = 9/10) +
  geom_point(data=pca[1, ], size=5,alpha = 9/10) +
  geom_point(data=pca[2, ], size=5,alpha = 9/10,)

b <- b + scale_x_reverse() + coord_equal()+ theme_classic() +scale_shape_manual(values=13:ka)
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
aca<-b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("pca_nean_chimp_shapes_new_axes_fix.png", aca, width=7, height=5, units="in", dpi=500)
getwd()



