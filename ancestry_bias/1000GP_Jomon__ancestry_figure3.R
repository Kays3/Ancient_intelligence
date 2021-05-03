#Ancestral bias in count polygenic score for 1000G people
#ALL snps used in P-link 
#1000GP 
#may 3, 2021
#Author kd

#plotting
library(dplyr)

library(ggplot2)

library(tidyverse)
library(ggpubr)
library(rstatix)





#2d plot for negative and positive 1000 Genome people
setwd("~/Jomon_intelligence")

GP_all<-read.csv("true_intelligence_1000G_J_12k.csv")
which(GP_all$subject=="Jomon")
jomon<-GP_all[1657,]

GP_all<-GP_all[-1657,]

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

GP_all_2 <- as_tibble(data.frame(GP_all, spp))
GP_all_2$spp <- factor(GP_all_2$spp)


levels(GP_all_2$spp)



effect1<- GP_all_2 %>% group_by(spp) %>% get_summary_stats(effect_sum,type = "mean_sd")
positive1<-GP_all_2 %>% group_by(spp) %>% get_summary_stats(positive_effect,type = "mean_sd")
negative1<-GP_all_2 %>% group_by(spp) %>% get_summary_stats(negative_effect,type = "mean_sd")

GP_all_3<-GP_all_2 %>%group_by(spp) %>% select(spp,effect_sum)

#no jomon just 1000 Genome people preview

ggboxplot(GP_all_2, x = "spp", y = "effect_sum")


dat<-GP_all_2
# ANOVA analysis
x <- which(names(dat) == "spp") # name of grouping variable
y <- which(
  names(dat) == "effect_sum" # names of variables to test
)
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("Europe", "Africa"),
                       c("Europe", "East_Asia"),
                       c("Europe", "Ad_Mixed_American"),
                       c("Europe", "South_Asia"),
                       c("East_Asia", "Ad_Mixed_American"),
                       c("Africa", "Ad_Mixed_American"),
                       c("South_Asia", "Africa"),
                       c("South_Asia", "East_Asia")
) 
# comparisons for post-hoc tests
# Edit until here
# Edit at your own risk
#library(ggpubr)
for (i in y) {
  for (j in x) {
    p <- ggboxplot(dat,
                   x = colnames(dat[j]), y = colnames(dat[i]),
                   color = colnames(dat[j]),
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
    )
    print(
      p + stat_compare_means(aes(label = paste0(..method.., ", p-value = ", ..p.format..)),
                             method = method1, label.y = max(dat[, i], na.rm = TRUE)
      )
      + stat_compare_means(comparisons = my_comparisons, method = method2, label = "p.format") # remove if p-value of ANOVA or Kruskal-Wallis test >= alpha
    )
  }
}



# Compute the analysis of variance
res.aov <- aov(effect_sum ~ spp, data = GP_all_2)
# Summary of the analysis
summary(res.aov)

#post-hoc testing in detail

lal<-TukeyHSD(res.aov,)
bogo<-as.data.frame(lal$spp)
write.csv(bogo,"tukey_12k_1000gp.csv",row.names = TRUE)

pairwise.t.test(GP, ses, p.adj = "bonf")


