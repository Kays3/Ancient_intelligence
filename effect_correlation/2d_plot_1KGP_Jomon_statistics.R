#absolute count of alleles plot for 1000G people
#ALL snps used in P-link 
#1000GP only
#oct 22, 2020
#Author kd

#plotting
rm(list = ls())
library("ggExtra")
library("ggplot2")
library(dplyr)


setwd("~/effect_correlation")

#2d plot for neg and positive ASIA

GP_all<-read.csv("merged_1G.csv")

GP2<-data.frame(rep('GP_all',length(GP_all$subject)),GP_all$positive_count,GP_all$negative_count)
colnames(GP2)<-c("subject","positive_count","negative_count")
head(GP2)
#high counts
GP3<-GP_all[GP_all$total_ALL_1000GP>240,]
#low counts
GP5<-GP_all[GP_all$total_ALL_1000GP<196,]

#top negative 
GP6<-GP_all[GP_all$negative_count>132,]


ggplot(GP3, aes(x = factor(subject))) +
  geom_bar()


GP4<-as.data.frame(table(GP3$subject))

GP4 %>%
  mutate(Var1=factor(Var1, Var1)) %>%
  ggplot( aes(x=Var1, y=Freq) ) +
  geom_segment( aes(x=Var1 ,xend=Var1, y=0, yend=Freq), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme_ipsum(base_family = "Times New Roman",
              base_size = 15,) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("")

#ggsave("Figure_3C.png", width = 30, height = 17, units = "cm",dpi = "retina")



#Geogrpahical highlights

# remake data.frame
spp <- rep(NA, length(GP_all$subject))

spp[grep("CHB", GP_all$subject)] <- "East_Asia"
spp[grep("JPT", GP_all$subject)] <- "East_Asia"
spp[grep("CHS", GP_all$subject)] <- "East_Asia"
spp[grep("CDX", GP_all$subject)] <- "East_Asia"
spp[grep("KHV", GP_all$subject)] <- "East_Asia"
spp[grep("CEU", GP_all$subject)] <- "Europe"
spp[grep("TSI", GP_all$subject)] <- "Europe"
spp[grep("FIN", GP_all$subject)] <- "Europe"
spp[grep("GBR", GP_all$subject)] <- "Europe"
spp[grep("IBS", GP_all$subject)] <- "Europe"
spp[grep("YRI", GP_all$subject)] <- "Africa"
spp[grep("LWK", GP_all$subject)] <- "Africa"
spp[grep("GWD", GP_all$subject)] <- "Africa"
spp[grep("MSL", GP_all$subject)] <- "Africa"
spp[grep("ESN", GP_all$subject)] <- "Africa"
spp[grep("ASW", GP_all$subject)] <- "Africa"
spp[grep("ACB", GP_all$subject)] <- "Africa"
spp[grep("MXL", GP_all$subject)] <- "Ad_Mixed_American"
spp[grep("PUR", GP_all$subject)] <- "Ad_Mixed_American"
spp[grep("CLM", GP_all$subject)] <- "Ad_Mixed_American"
spp[grep("PEL", GP_all$subject)] <- "Ad_Mixed_American"
spp[grep("GIH", GP_all$subject)] <- "South_Asia"
spp[grep("PJL", GP_all$subject)] <- "South_Asia"
spp[grep("BEB", GP_all$subject)] <- "South_Asia"
spp[grep("STU", GP_all$subject)] <- "South_Asia"
spp[grep("ITU", GP_all$subject)] <- "South_Asia"
spp[grep("AltaiNea", GP_all$subject)] <- "AltaiNea"
spp[grep("Denisovan", GP_all$subject)] <- "Denisovan"
spp[grep("Jomon", GP_all$subject)] <- "Jomon"
spp[grep("taweh", GP_all$subject)] <- "Chimpanzee"



new <- as_tibble(data.frame(GP_all, spp))
tail(new)
new$spp <- factor(new$spp)

colnames(new)

EU<-new[new$spp =="Europe",]
AF<-new[new$spp =="Africa",]


boblo<-rbind(EU,AF)
cor.test(boblo$positive_count,boblo$negative_count)
#tri
#correlation positive and negative GEO
b <- ggplot(boblo, aes(negative_count,positive_count, col = spp)) +
  geom_point(size = 5,alpha = 5/10)


b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
#ggsave("2D_GEO_all.png", aca, width=7, height=5, units="in", dpi=500)


#correlation positive and negative GEO
b <- ggplot(new, aes(negative_count,positive_count, col = spp)) +
  geom_point(size = 5,alpha = 5/10)

  
b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
#ggsave("2D_GEO_all.png", aca, width=7, height=5, units="in", dpi=500)


#correlation positive and negative EU
b <- ggplot(EU, aes(negative_count,positive_count, col = spp)) +
  geom_point(size = 5,alpha = 5/10)


b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
#ggsave("2D_GEO_all.png", aca, width=7, height=5, units="in", dpi=500)


lm_fit_EU <- lm(negative_count ~ positive_count, data=EU)
summary(lm_fit_EU)

#label top 1 % highest and top 1 % lowest
# filter dataframe to get data to be highligheted
highlight_high <- new %>% 
  filter(total_ALL_1000GP>=240)


highlight_low <- new %>% 
  filter(total_ALL_1000GP<=196)

b <- ggplot(new, aes(negative_count,positive_count)) +
  geom_point(size = 3,alpha = 5/10)+
  geom_point(data=highlight_high, aes(x=negative_count,y=positive_count), color='red',size=3)+
  geom_point(data=highlight_low, aes(x=negative_count,y=positive_count), color='blue',size=3)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1), colour="black",se = FALSE)

b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))

#ggsave("2D_high_low_margins_1_percent.png", aca, width=7, height=5, units="in", dpi=500)


#try eu AFR
b <- ggplot(boblo, aes(negative_count,positive_count,col = spp)) +
  geom_point(size = 3,alpha = 5/10)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1), colour="black",se = FALSE)

b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))

ggsave("2D_high_positive_cor.png", aca, width=7, height=5, units="in", dpi=500)
getwd()

#all GEO new
b <- ggplot(new, aes(negative_count,positive_count,col = spp)) +
  geom_point(size = 3,alpha = 5/10)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1), colour="black",se = FALSE)

b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))

ggsave("2D_all_GEO_cor.png", aca, width=7, height=5, units="in", dpi=500)
getwd()
#plot 1 % red


b <- ggplot(highlight_high, aes(negative_count,positive_count)) +
  geom_point(size = 3,alpha = 10/10, color = 'red')+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1), colour="black",se = FALSE)


b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
ggsave("2D_high_top1_percent.png", aca, width=7, height=5, units="in", dpi=500)

#plot top 1 % GEO
#highlight Europe

non_c <- setdiff(levels(highlight_high$spp), "Europe")
levels(highlight_high$spp) <- list(Europe = "Europe", "World" = non_c)



qu1<-ggplot(highlight_high, aes(x=as.factor(spp), fill=as.factor(spp) )) + 
  geom_bar( ) +
  scale_fill_manual(values = c("red","red") ) + 
  theme_classic()+
  labs(
    x = "Top 1 % scores",
    y = "Count")
qu1
  ggsave("top1_percent_GEO_EU.png", qu1, width=7, height=5, units="in", dpi=500)



b <- ggplot(highlight_high, aes(negative_count,positive_count, col = spp)) +
  geom_point(size = 3,alpha = 10/10)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black",se = FALSE)
  
  
b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
ggsave("2D_high_top1_percent_GEO_EU.png", aca, width=7, height=5, units="in", dpi=500)

lm_fit_high <- lm(negative_count ~ positive_count, data=highlight_high)
summary(lm_fit_high)

lm_fit_low <- lm(negative_count ~ positive_count, data=highlight_low)
summary(lm_fit_low)
cor.test(highlight_high$negative_count,highlight_high$positive_count,method="pearson")
cor.test(highlight_low$negative_count,highlight_low$positive_count,method="pearson")

#plot bottom 1 %
#plot bottom 1 % blue

b <- ggplot(highlight_low, aes(negative_count,positive_count)) +
  geom_point(size = 3,alpha = 10/10, color = 'blue')+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1), colour="black",se = FALSE)


b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
ggsave("2D_high_bottom1_percent.png", aca, width=7, height=5, units="in", dpi=500)

#GEO


qu2<-ggplot(highlight_low, aes(x=as.factor(spp), fill=as.factor(spp) )) + 
  geom_bar( ) +
  scale_fill_manual(values = c("blue","blue") ) + 
  theme(legend.position="none")+
  labs(
    x = "Bottom 1 % scores",
    y = "Count")
qu2

ggsave("bottom1_percent_GEO.png", qu2, width=7, height=5, units="in", dpi=500)
 
b <- ggplot(highlight_low, aes(negative_count,positive_count, col = spp)) +
  geom_point(size = 3,alpha = 10/10)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,
              aes(group=1),colour="black",se = FALSE)


b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Positive alleles"))
ggsave("2D_high_bottom1_percent_GEO1.png", aca, width=7, height=5, units="in", dpi=500)

lm_fit_low <- lm(negative_count ~ positive_count, data=highlight_low)
summary(lm_fit_low)

#ggsave("corr1", aca, width=7, height=5, units="in", dpi=500)

#correlation positive and total_GEO
b <- ggplot(new, aes(positive_count, total_ALL_1000GP, col = spp)) +
  geom_point(size = 5,alpha = 5/10)

b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Positive alleles")) + ylab(paste0("Total alleles "))

aca<-b + xlab(paste0("Positive alleles")) + ylab(paste0("Total alleles "))


#correlation negative and total
b <- ggplot(new, aes(negative_count, total_ALL_1000GP, col = spp)) +
  geom_point(size = 5,alpha = 5/10)

b <- b + coord_equal()+ theme_classic() 
b + xlab(paste0("Negative alleles")) + ylab(paste0("Total alleles "))

aca<-b + xlab(paste0("Negative alleles")) + ylab(paste0("Total alleles "))



df1 <- data.frame(x = GP2$positive_count, y = GP2$negative_count)
colnames(df1)<-c("positive", "negative")

p1 <- ggplot(df1, aes(positive, negative)) + geom_point() + theme_bw()+
  geom_smooth(method = "lm", se = FALSE)
p1
ggMarginal(p1, size = 2 , type = "histogram",
           col = "black", fill = "gray",
)

cor.test(boblo$positive_count,boblo$negative_count)
