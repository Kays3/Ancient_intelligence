# load tidyverse package
#add Jomon
#1402 snsps work  - NO GOMI Afanasievo included
#june 2 2021  - style for paper 12 k snps total updated
# 1402 out of 1880 with weird pca -- cloudyu

rm(list = ls())
library(tidyverse)
#library(wesanderson)
#library(viridis)

library("ggsci")
setwd("~/testdir/savage/plink_intelligence/Afanasievo/mod_1000g_jom_afa/")

library(plyr)
library(readr)

#load 1000gp phase 3 list of people
phase3<- read.table("~/testdir/savage/phase3_people.csv", 
                    sep = "\t" , header = F,
                    na.strings ="", stringsAsFactors= F, 
)

#merged PCA
#1000G_jomon_1402_AFA


# read in data
pca <- read_table2("./merged_update_1000G_J_Afa_1880.eigenvec", col_names = FALSE)
eigenval <- scan("./merged_update_1000G_J_Afa_1880.eigenval")



# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
names(pca)[1] <- "ind"
# set names per phase3


names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
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
m1 <- merge(pca, pops, by.x = "ind",all.x = TRUE)
pca<-m1

tail(pca)

#find Jomon
which(pca$popa == "Jomon")
#1661
which(pca$popa == "Afanasievo_s1")
#1658
which(pca$popa == "Afanasievo_s2")
#1660
which(pca$popa == "Afanasievo_m")
#1657
which(pca$popa == "Afanasievo_f")
#1659

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
spp[grep("Afanasievo_s2", pca$popa)] <- "Afanasievo"
spp[grep("Afanasievo_s1", pca$popa)] <- "Afanasievo"
spp[grep("Afanasievo_m", pca$popa)] <- "Afanasievo"
spp[grep("Afanasievo_f", pca$popa)] <- "Afanasievo"

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



#new pca
setwd("/Users/kaisar_dauyey/testdir/savage/plink_intelligence/results/")


ka<- c(15+ nlevels(pca$spp))

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = spp )) +
  geom_point(size = 3,alpha = 7/10) +
  geom_point(data=pca[1661, ], size=9,alpha = 10/10) +
  geom_point(data=pca[1658, ], size=6,alpha = 9/10) +
  geom_point(data=pca[1657, ], size=6,alpha = 9/10) +
  geom_point(data=pca[1660, ], size=6,alpha = 9/10) +
  geom_point(data=pca[1659, ], size=6,alpha = 9/10) +
  scale_color_npg()  

b <- b + coord_equal()+ theme_classic() +scale_shape_manual(values=15:ka)
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
aca<-b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("1402_pca_1000G_jomon_Afa.png", aca, width=7, height=5, units="in", dpi=500)

getwd()


# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp)) +
  geom_point(size = 2,alpha = 7/10,shape =17) +
  geom_point(data=pca[1661, ], size=6,alpha = 9/10,shape =16) +
  geom_point(data=pca[1658, ], size=6,alpha = 9/10) +
  geom_point(data=pca[1657, ], size=6,alpha = 9/10) +
  geom_point(data=pca[1660, ], size=6,alpha = 9/10) +
  geom_point(data=pca[1659, ], size=6,alpha = 9/10) +
  scale_color_npg() 


b <- b + coord_equal()+ theme_classic()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
aca<-b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

ggsave("pca_1402_1000G_jomon_Afanasievo.png", aca, width=7, height=5, units="in", dpi=500)



#enough for paper....

