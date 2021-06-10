#june 2, 2021
#1402 snps common 1000G and Jomon _Afanasievo
#genotypes extracted
#polygenic score calculated
#EA proxy replication not for all snps...

#Author kaisar dauyey

rm(list = ls())
library(dplyr)

library(ggplot2)
library(cowplot)

library(tidyverse)
library(ggpubr)
library(rstatix)


setwd("/Users/kaisar_dauyey/testdir/savage/plink_intelligence/results")


#no gomi
GP_nogomi_j<-read_csv("~/testdir/savage/plink_intelligence/results/PGS_1000G_J_Afa_1402.csv")

GP3<- GP_nogomi_j %>% select(subject,effect_sum) %>%
  group_by(subject)
GP3<-GP3[-c(1657,1658,1659,1660,1661),]


jomon_afa<-GP_nogomi_j[c(1657,1658,1659,1660,1661),]




head(jomon_afa)



#####
#9k SNPs pipeline 
#Afanasievo

#POPULATION PARAMETER CALCULATIONS
pop_sd242 <- sd(GP3$effect_sum)*sqrt((length(GP3$effect_sum)-1)/(length(GP3$effect_sum)))
pop_mean242 <- mean(GP3$effect_sum)
# sd function in R uses the sample standard deviation
# and not the population standard deviation, difference for 


#z score estimation for F23
j242<-jomon_afa$effect_sum[5]
z242 <- (j242 - pop_mean242) / pop_sd242
p_yellow_j242 <- pnorm(z242,lower.tail = TRUE, log.p = FALSE)  



#z score estimation for Afa_mother
afa_m<-jomon_afa$effect_sum[1]
z_afa_m <- (afa_m - pop_mean242) / pop_sd242
p_yellow_afa_m <- pnorm(z_afa_m,lower.tail = TRUE, log.p = FALSE) 



#z score estimation for Afa_son1
afa_s1<-jomon_afa$effect_sum[2]
z_afa_s1 <- (afa_s1 - pop_mean242) / pop_sd242
p_yellow_afa_s1 <- pnorm(z_afa_s1,lower.tail = TRUE, log.p = FALSE) 


#z score estimation for Afa_father
afa_f<-jomon_afa$effect_sum[3]
z_afa_f <- (afa_f - pop_mean242) / pop_sd242
p_yellow_afa_f <- pnorm(z_afa_f,lower.tail = TRUE, log.p = FALSE) 

#z score estimation for Afa_son2
afa_s2<-jomon_afa$effect_sum[4]
z_afa_s2 <- (afa_s2 - pop_mean242) / pop_sd242
p_yellow_afa_s2 <- pnorm(z_afa_s2,lower.tail = TRUE, log.p = FALSE) 


#change labels per number appropriately
# ggplot  NORMAL distribution for 1000 GP vs F23 non-standrd
b<- ggplot(GP3, aes(x = GP3$effect_sum),) +
  stat_function(
    fun = dnorm,
    args = with(GP3, c(mean = pop_mean242, sd = pop_sd242))
  )+  
  geom_vline(aes(xintercept=jomon_afa$effect_sum[5],
                 color="Funadomari Jomon"),     linetype="solid", size=3)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[1],
                 color="Afanasievo Mother"),     linetype="solid", size=3)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[2],
                 color="Afanasievo Son1"),     linetype="solid", size=3)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[3],
                 color="Afanasievo Father"),     linetype="solid", size=3)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[4],
                 color="Afanasievo Son2"),     linetype="solid", size=3)+
  
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
  
  
  scale_color_manual(name = "Polygenic score, 1402 SNPs", values = c("Funadomari Jomon" = "#33FFFF",
                                                                     "Afanasievo Mother" = "red",
                                                                     "Afanasievo Father" = "blue",
                                                                     "Afanasievo Son1" = "green",
                                                                     "Afanasievo Son2" = "purple",
                                                                     "Mean 1000 Genomes Project" = "black"))
b
#ggsave("1000GP_jomon_afa_9128_2.pdf",b,width = 30, height = 15, units = "cm", dpi = "retina")
ggsave("1000GP_jomon_afa_1402.png",b, width = 30, height = 15, units = "cm",dpi = "retina")

#####
