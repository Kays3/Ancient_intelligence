#june 11, 2021
#1402 snps common 1000G and Jomon _Afanasievo
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
library(extrafont)
#font_import()
#loadfonts()

#setwd("/data/")
setwd("/Users/kaisar_dauyey/testdir/savage/plink_intelligence/scripts/paper/data")

#no gomi
GP_nogomi_j<-read_csv("PGS_1000G_J_Afa_1402.csv")

GP3<- GP_nogomi_j %>% select(subject,effect_sum) %>%
  group_by(subject)
GP3<-GP3[-c(1657,1658,1659,1660,1661),]


jomon_afa<-GP_nogomi_j[c(1657,1658,1659,1660,1661),]




head(jomon_afa)



#####
#1402 SNPs pipeline 
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
                 color="Funadomari Jomon"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[1],
                 color="Afanasievo Mother"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[2],
                 color="Afanasievo Son1"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[3],
                 color="Afanasievo Father"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa$effect_sum[4],
                 color="Afanasievo Son2"),     linetype="solid", size=1)+
  
  geom_vline(aes(xintercept=median(GP3$effect_sum),
                 color="Mean 1000 Genomes Project"), linetype="solid", size=1)+
  #Set the size of the plotting window
  theme_pubclean()+
  
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 7, vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 7, vjust = 1),
    plot.caption = element_text(angle = 0, size = 7, vjust = 1),
    
    axis.text.x = element_text(angle = 45, size = 7, hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 7, vjust = 0.5),
    axis.title = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    
    #Legend
    legend.key = element_blank(),       #removes the border
    legend.key.size = unit(0.3, 'cm'),        #Sets overall area/size of the legend
    legend.text = element_text(size = 5),   #Text size
    title=element_text(size = 5),
  )+
  #Change the size of the icons/symbols in the legend
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  #Set x- and y-axes labels
  xlab('Score') +
  ylab('Density') +
  
  
  scale_color_manual(name = "Polygenic score, 1402 SNPs", values = c("Funadomari Jomon" = "#4DBBD5FF",
                                                                     "Afanasievo Mother" = "#E64B35FF",
                                                                     "Afanasievo Father" = "#3C5488FF",
                                                                     "Afanasievo Son1" = "#00A087FF",
                                                                     "Afanasievo Son2" = "#8491B499",
                                                                     "Mean 1000 Genomes Project" = "#7E6148FF"))
b
# ggsave(filename = "~/Desktop/Figure1A.pdf", plot = b, width = 17.5, height = 10, dpi = 300, units = "cm")

#9k work

#data loading
GP_9k<-read_csv("PGS_1000G_J_Afa_9k.csv")

GP3_9k<- GP_9k %>% select(subject,effect_sum) %>%
  group_by(subject)
GP3_9k<-GP3_9k[-c(1657,1658,1659,1660,1661),]


jomon_afa_9k<-GP_9k[c(1657,1658,1659,1660,1661),]




head(jomon_afa_9k)



#####
#9k SNPs pipeline 
#Afanasievo

#POPULATION PARAMETER CALCULATIONS
pop_sd_9k <- sd(GP3_9k$effect_sum)*sqrt((length(GP3_9k$effect_sum)-1)/(length(GP3_9k$effect_sum)))
pop_mean_9k <- mean(GP3_9k$effect_sum)
# sd function in R uses the sample standard deviation
# and not the population standard deviation, difference for 


#z score estimation for F23
j9k<-jomon_afa_9k$effect_sum[5]
z9k <- (j9k - pop_mean_9k) / pop_sd_9k
p_yellow_j9k <- pnorm(z9k,lower.tail = TRUE, log.p = FALSE)  



#z score estimation for Afa_mother
afa_m_9k<-jomon_afa_9k$effect_sum[1]
z_afa_m_9k <- (afa_m_9k - pop_mean_9k) / pop_sd_9k
p_yellow_afa_m_9k <- pnorm(z_afa_m_9k,lower.tail = TRUE, log.p = FALSE) 



#z score estimation for Afa_son1
afa_s1_9k<-jomon_afa_9k$effect_sum[2]
z_afa_s1_9k <- (afa_s1_9k - pop_mean_9k) / pop_sd_9k
p_yellow_afa_s1_9k <- pnorm(z_afa_s1_9k,lower.tail = TRUE, log.p = FALSE) 


#z score estimation for Afa_father
afa_f_9k<-jomon_afa_9k$effect_sum[3]
z_afa_f_9k <- (afa_f_9k - pop_mean_9k) / pop_sd_9k
p_yellow_afa_f_9k <- pnorm(z_afa_f_9k,lower.tail = TRUE, log.p = FALSE) 

#z score estimation for Afa_son2
afa_s2_9k<-jomon_afa_9k$effect_sum[4]
z_afa_s2_9k <- (afa_s2_9k - pop_mean_9k) / pop_sd_9k
p_yellow_afa_s2_9k <- pnorm(z_afa_s2_9k,lower.tail = TRUE, log.p = FALSE) 



# ggplot  NORMAL distribution for 1000 GP vs F23 non-standrd
b2<- ggplot(GP3_9k, aes(x = GP3_9k$effect_sum),) +
  stat_function(
    fun = dnorm,
    args = with(GP3_9k, c(mean = pop_mean_9k, sd = pop_sd_9k))
  )+  
  geom_vline(aes(xintercept=jomon_afa_9k$effect_sum[5],
                 color="Funadomari Jomon"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa_9k$effect_sum[1],
                 color="Afanasievo Mother"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa_9k$effect_sum[2],
                 color="Afanasievo Son1"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa_9k$effect_sum[3],
                 color="Afanasievo Father"),     linetype="solid", size=1)+
  geom_vline(aes(xintercept=jomon_afa_9k$effect_sum[4],
                 color="Afanasievo Son2"),     linetype="solid", size=1)+
  
  geom_vline(aes(xintercept=median(GP3_9k$effect_sum),
                 color="Mean 1000 Genomes Project"), linetype="solid", size=1)+
  #Set the size of the plotting window
  theme_pubclean()+
  
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size =  7,   vjust = 1),
    plot.subtitle = element_text(angle = 0, size =  7,   vjust = 1),
    plot.caption = element_text(angle = 0, size =  7,   vjust = 1),
    
    axis.text.x = element_text(angle = 45, size =  7,   hjust = 1.10),
    axis.text.y = element_text(angle = 0, size =  7,   vjust = 0.5),
    axis.title = element_text(size =  7  ),
    axis.title.x = element_text(size =  7 ),
    axis.title.y = element_text(size =  7  ),
    
    #Legend
    legend.key = element_blank(),       #removes the border
    legend.key.size = unit(0.3, 'cm'),        #Sets overall area/size of the legend
    legend.text = element_text(size = 5),   #Text size
    title=element_text(size = 5),
  )+
  #Change the size of the icons/symbols in the legend
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  #Set x- and y-axes labels
  xlab('Score') +
  ylab('Density') +
  
  
  scale_color_manual(name = "Polygenic score, 9128 SNPs", values = c("Funadomari Jomon" = "#4DBBD5FF",
                                                                     "Afanasievo Mother" = "#E64B35FF",
                                                                     "Afanasievo Father" = "#3C5488FF",
                                                                     "Afanasievo Son1" = "#00A087FF",
                                                                     "Afanasievo Son2" = "#8491B499",
                                                                     "Mean 1000 Genomes Project" = "#7E6148FF"))

b2

Figure1 <- plot_grid(b, b2, labels = "AUTO", ncol = 1, align = 'v')

Figure1

ggsave("~/Desktop/Figure1.pdf", Figure1, width=17.5, height=15, units="cm", dpi=300)

