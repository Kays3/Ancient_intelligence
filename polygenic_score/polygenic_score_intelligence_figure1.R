#May 3 , 2021
#12 037 common in Jomon and 1000 Genome subjects were used
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


setwd("~/Jomon_intelligence")

#loading data
GP_nogomi_j<-read_csv("true_intelligence_1000G_J_12k.csv")

GP3<- GP_nogomi_j %>% select(subject,effect_sum) %>%
  group_by(subject)
GP3<-GP3[-c(1657),]

#subset jomon from 1000 genomes
jomon242<-GP_nogomi_j[1657,]

#12k SNPs pipeline 


#POPULATION PARAMETER CALCULATIONS
pop_sd242 <- sd(GP3$effect_sum)*sqrt((length(GP3$effect_sum)-1)/(length(GP3$effect_sum)))
pop_mean242 <- mean(GP3$effect_sum)
# sd function in R uses the sample standard deviation
# and not the population standard deviation, difference for 


#z score estimation for F23_242
j242<-jomon242$effect_sum[1]
z242 <- (j242 - pop_mean242) / pop_sd242
p_yellow_j242 <- pnorm(z242,lower.tail = TRUE, log.p = FALSE)  




# ggplot  NORMAL distribution for 1000 GP vs F23 non-standrd
b<- ggplot(GP3, aes(x = GP3$effect_sum),) +
  stat_function(
    fun = dnorm,
    args = with(GP3, c(mean = pop_mean242, sd = pop_sd242))
  )+  
  geom_vline(aes(xintercept=jomon242$effect_sum[1],
                 color="Funadomari Jomon"),     linetype="dashed", size=3)+
  
  geom_vline(aes(xintercept=median(GP3$effect_sum),
                 color="Mean 1000 Genomes Project"), linetype="solid", size=3)+
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
  
  
  scale_color_manual(name = "Polygenic score, 12037 SNPs", values = c("Funadomari Jomon" = "#33FFFF",
                                                                      "Mean 1000 Genomes Project" = "black"))
#if you want to safe plots

ggsave("1000GP_jomon_12037.pdf",width = 30, height = 15, units = "cm", dpi = "retina")
ggsave("1000GP_jomon_12037.png", width = 30, height = 15, units = "cm",dpi = "retina")

#####
