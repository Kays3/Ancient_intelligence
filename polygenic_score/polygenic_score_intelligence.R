#September 30, 2020
#72 snps common in Neanderthal and Denisovan were compared to 1000G and Jomon and CHIMP
#242 common in Jomon and 1000 Genome subjects were used as well
#genotypes extracted
#polygenic score calculated from GWAS meta-analysis by Savage et al 2018
#normalization of scores omitted for the purpose of data comparison between groups
#Author kaisar dauyey

rm(list = ls())
library(dplyr)

library(ggplot2)
library(cowplot)

library(tidyverse)
library(ggpubr)
library(rstatix)


setwd("~/polygenic_score")

nean<-read.csv("intelligence_nean_deni.csv")
#72
GP_J_chimp_nean_deni<-read_csv("intelligence_1000G_J_chimp_nean_deni72.csv")

#242
GP_chimp_j<-read_csv("intelligence_chimp.csv")

GP3<- GP_chimp_j %>% select(subject,effect_sum) %>%
  group_by(subject)
GP3<-GP3[-c(1657,2506),]

buddies_242<-rbind(GP_chimp_j[1657,], GP_chimp_j[2506,])
jomon242<-GP_chimp_j[1657,]


#72 SNPs pipeline
head(nean)

GP2<- GP_J_chimp_nean_deni %>% select(subject,effect_sum) %>%
  group_by(subject)
GP2a<-GP2[-c(1,2,1659,2508),]
buddies<-GP2[c(1,2,1659,2508),]
ggboxplot(GP2, x = "subject", y = "effect_sum")


head(buddies)



ggboxplot(GP2, x = "subject", y = "effect_sum")+
  geom_hline(aes(yintercept=buddies$effect_sum[1],
                 color="Jomon_F23"), linetype="dashed",
             size=1)+
  geom_hline(aes(yintercept=buddies$effect_sum[2],
                 color="Neanderthal"), linetype="dashed",
             size=1)+
  geom_hline(aes(yintercept=buddies$effect_sum[3],
                 color="Denisovan"), linetype="dashed",
             size=1)+
  geom_hline(aes(yintercept=buddies$effect_sum[4],
                 color="Chimpanzee"), linetype="dashed",
             size=1)+
  geom_hline(aes(yintercept=median(GP2$effect_sum),
                 color="Median"), linetype="dashed", size=1)

ggsave("1000GP_populations_merged_nean_chimp_prelim2.png", width = 30, height = 13, units = "cm",dpi = "retina")


#POPULATION PARAMETER CALCULATIONS
pop_sd <- sd(GP2a$effect_sum)*sqrt((length(GP2a$effect_sum)-1)/(length(GP2$effect_sum)))
pop_mean <- mean(GP2a$effect_sum)
# sd function in R uses the sample standard deviation
# and not the population standard deviation, difference for 
#z score estimation for F23
j1<-buddies$effect_sum[1]
z1 <- (j1 - pop_mean) / pop_sd
p_yellow1 <- 2*pnorm(z1)  

p_yellow1<-pnorm(z1 , lower.tail = TRUE, log.p = FALSE)

#nean
j2<-buddies$effect_sum[2]
z2 <- (j2 - pop_mean) / pop_sd
p_yellow2 <- pnorm(z2,lower.tail = TRUE, log.p = FALSE)  

#deni
j3<-buddies$effect_sum[3]
z3 <- (j3 - pop_mean) / pop_sd
p_yellow3 <- pnorm(z3, lower.tail = TRUE, log.p = FALSE)  

#chimp
j4<-buddies$effect_sum[4]
z4 <- (j4 - pop_mean) / pop_sd
p_yellow4 <- pnorm(z4,lower.tail = FALSE, log.p = FALSE)  

#bonferroni adj p values  - no need
#pibi<-c(p_yellow1,p_yellow2,p_yellow3,p_yellow4)
#pibi_ad<-round(p.adjust(pibi, "bonferroni" ),3)
#pibi_ad

# ggplot  NORMAL distribution for 1000 GP vs F23_72_nean_deni_CHIMP
ggplot(GP2a, aes(x = GP2a$effect_sum),) +
  stat_function(
    fun = dnorm,
    args = with(GP2a, c(mean = pop_mean, sd = pop_sd))
  )+  
  geom_vline(aes(xintercept=buddies$effect_sum[1],
                 color="Funadomari Jomon"),     linetype="dashed", size=3)+
  geom_vline(aes(xintercept=buddies$effect_sum[2],
                 color="Neanderthal"),     linetype="dashed", size=3)+
  geom_vline(aes(xintercept=buddies$effect_sum[3],
                 color="Denisovan"),     linetype="dashed", size=3)+
  geom_vline(aes(xintercept=buddies$effect_sum[4],
                 color="Chimpanzee"),     linetype="dashed", size=3)+
  geom_vline(aes(xintercept=mean(GP2$effect_sum),
                 color="Mean 1000 Genome Project"), linetype="dashed", size=3)+
  #Set the size of the plotting window
  theme_pubclean()+
  
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 45, size = 16, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 16, face = 'bold'),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    
    #Legend
    legend.key = element_blank(),       #removes the border
    legend.key.size = unit(1, 'cm'),        #Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = 'bold'),   #Text size
    title=element_text(size = 14, face = 'bold'),
  )+
  #Change the size of the icons/symbols in the legend
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  #Set x- and y-axes labels
  xlab('Score') +
  ylab('Density') +
  
  
  scale_color_manual(name = "Individual score, 72 SNPs", values = c("Funadomari Jomon" = "#33FFFF",
                                                                    "Neanderthal" = "#FFCC00",
                                                           "Denisovan" = "#66CC00",
                                                           "Chimpanzee" = "#333333",
                                                           "Mean 1000 Genome Project" = "#CA1036"))

ggsave("1000GP_72_new_chimp.pdf",width = 30, height = 15, units = "cm", dpi = "retina")
ggsave("1000GP_72_new_chimp.png", width = 30, height = 15, units = "cm",dpi = "retina")

#####
#242 SNPs pipeline


#POPULATION PARAMETER CALCULATIONS
pop_sd242 <- sd(GP3$effect_sum)*sqrt((length(GP3$effect_sum)-1)/(length(GP3$effect_sum)))
pop_mean242 <- mean(GP3$effect_sum)
# sd function in R uses the sample standard deviation
# and not the population standard deviation, difference for 


#z score estimation for F23_242
j242<-buddies_242$effect_sum[1]
z242 <- (j242 - pop_mean242) / pop_sd242
p_yellow_j242 <- pnorm(z242,lower.tail = TRUE, log.p = FALSE)  


chimp242<-buddies_242$effect_sum[2]
z_chimp242 <- (chimp242 - pop_mean242) / pop_sd242
p_yellow_chimp242 <- pnorm(z_chimp242,lower.tail = FALSE, log.p = FALSE) 
# comparison no adjustment



# ggplot  NORMAL distribution for 1000 GP vs F23 non-standrd
ggplot(GP3, aes(x = GP3$effect_sum),) +
  stat_function(
    fun = dnorm,
    args = with(GP3, c(mean = pop_mean242, sd = pop_sd242))
  )+  
  geom_vline(aes(xintercept=buddies_242$effect_sum[1],
                 color="Funadomari Jomon"),     linetype="dashed", size=3)+
  geom_vline(aes(xintercept=buddies_242$effect_sum[2],
                 color="Chimpanzee"),     linetype="dashed", size=3)+
  geom_vline(aes(xintercept=median(GP3$effect_sum),
                 color="Mean 1000 Genome Project"), linetype="dashed", size=3)+
  #Set the size of the plotting window
  theme_pubclean()+
  
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 45, size = 16, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 16, face = 'bold'),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    
    #Legend
    legend.key = element_blank(),       #removes the border
    legend.key.size = unit(1, 'cm'),        #Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = 'bold'),   #Text size
    title=element_text(size = 14, face = 'bold'),
  )+
  #Change the size of the icons/symbols in the legend
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  #Set x- and y-axes labels
  xlab('Score') +
  ylab('Density') +
  
 
  scale_color_manual(name = "Individual score, 242 SNPs", values = c("Funadomari Jomon" = "#33FFFF",
                                                                     "Chimpanzee" = "#333333",
                                                                     "Mean 1000 Genome Project" = "#CA1036"))

ggsave("1000GP_jomon_242_new2.pdf",width = 30, height = 15, units = "cm", dpi = "retina")
ggsave("1000GP_jomon_242_new2.png", width = 30, height = 15, units = "cm",dpi = "retina")

#####
