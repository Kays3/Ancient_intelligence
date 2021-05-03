#ALL snps used in P-link 
#1000GP 
#may 3, 2021
#PCA visualization
#Author kd
rm(list = ls())
library(tidyverse)


library(plyr)
library(readr)
setwd("~/Jomon_intelligence/pca")

#load 1000gp phase 3 list of people
phase3<- read.table("phase3_people.csv", 
                    sep = "\t" , header = F,
                    na.strings ="", stringsAsFactors= F, 
)

#merged PCA
#1000G_jomon_12k


# read in data
pca <- read_table2("./1000G_Jomon_no_gomi_2.eigenvec", col_names = FALSE)
eigenval <- scan("./1000G_Jomon_no_gomi_2.eigenval")



# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
names(pca)[1] <- "ind"
# set names per phase3


names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pops<-data.frame(phase3$V2,phase3$V3,stringsAsFactors=FALSE)
extra1<-c('Jomon', "Jomon")


pops<-rbind(pops,extra1)

tail(pops)
names(pops)<-c("ind","popa")
pops<-data.frame(pops)
m1 <- merge(pca, pops, by.x = "ind",all.x = TRUE)
pca<-m1

tail(pca)

#find Jomon
which(pca$popa == "Jomon")
#1657

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
spp[grep("Jomon", pca$popa)] <- "Jomon"
#spp[grep("taweh", pca$popa)] <- "Chimpanzee"

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
#to safe
#ggsave("pve4.png", ala, width=7, height=5, units="in", dpi=500)



#build pca plot


ka<- c(15+ nlevels(pca$spp))

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = spp )) +
  geom_point(size = 3,alpha = 1/10) +
  geom_point(data=pca[1657, ], size=5,alpha = 9/10)

b <- b + coord_equal()+ theme_classic() +scale_shape_manual(values=15:ka)
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
aca<-b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

#to safe edit comment
#ggsave("12k_pca_1000G_jomon_NO_GOMI.png", aca, width=7, height=5, units="in", dpi=500)



