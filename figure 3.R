#july 17, 2021
#no normality done since we need to do data comparison between groups
#12k racial differences in PGS
#modified
#Author Kaisar Dauyey
#uncomment to save plots


rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(pairwiseComparisons)
library(ggsignif)



#setwd("data/")
setwd("/Users/kaisar_dauyey/testdir/savage/plink_intelligence/scripts/paper/data")

#2d plot for neg and positive ASIA

GP_all<-read.csv("true_intelligence_1000G_J_12k.csv")
which(GP_all$subject=="Jomon")
jomon<-GP_all[1657,]

GP_all<-GP_all[-1657,]

spp <- rep(NA, length(GP_all$subject))

spp[grep("CHB", GP_all$subject)] <- "EAS"
spp[grep("JPT", GP_all$subject)] <- "EAS"
spp[grep("CHS", GP_all$subject)] <- "EAS"
spp[grep("CDX", GP_all$subject)] <- "EAS"
spp[grep("KHV", GP_all$subject)] <- "EAS"
spp[grep("CEU", GP_all$subject)] <- "EUR"
spp[grep("TSI", GP_all$subject)] <- "EUR"
spp[grep("FIN", GP_all$subject)] <- "EUR"
spp[grep("GBR", GP_all$subject)] <- "EUR"
spp[grep("IBS", GP_all$subject)] <- "EUR"
spp[grep("YRI", GP_all$subject)] <- "AFR"
spp[grep("LWK", GP_all$subject)] <- "AFR"
spp[grep("GWD", GP_all$subject)] <- "AFR"
spp[grep("MSL", GP_all$subject)] <- "AFR"
spp[grep("ESN", GP_all$subject)] <- "AFR"
spp[grep("ASW", GP_all$subject)] <- "AFR"
spp[grep("ACB", GP_all$subject)] <- "AFR"
spp[grep("MXL", GP_all$subject)] <- "AMR"
spp[grep("PUR", GP_all$subject)] <- "AMR"
spp[grep("CLM", GP_all$subject)] <- "AMR"
spp[grep("PEL", GP_all$subject)] <- "AMR"
spp[grep("GIH", GP_all$subject)] <- "SAS"
spp[grep("PJL", GP_all$subject)] <- "SAS"
spp[grep("BEB", GP_all$subject)] <- "SAS"
spp[grep("STU", GP_all$subject)] <- "SAS"
spp[grep("ITU", GP_all$subject)] <- "SAS"

GP_all_2 <- as_tibble(data.frame(GP_all, spp))
GP_all_2$spp <- factor(GP_all_2$spp)

#GP_all_2 %>% sample_n_by(spp, size = 1)
levels(GP_all_2$spp)
#GP_all_2$effect_sum<-scale(GP_all_2$effect_sum)



effect1<- GP_all_2 %>% group_by(spp) %>% get_summary_stats(effect_sum,type = "mean_sd")
positive1<-GP_all_2 %>% group_by(spp) %>% get_summary_stats(positive_effect,type = "mean_sd")
negative1<-GP_all_2 %>% group_by(spp) %>% get_summary_stats(negative_effect,type = "mean_sd")

GP_all_3<-GP_all_2 %>%group_by(spp) %>% select(spp,effect_sum)

#no jomon compariosn

ggboxplot(GP_all_2, x = "spp", y = "effect_sum")
#ggsave("1000GP_populations_12k_no_F23.png", width = 30, height = 13, units = "cm",dpi = "retina",)


dat<-GP_all_2
# Edit from here
x <- which(names(dat) == "spp") # name of grouping variable
y <- which(
  names(dat) == "effect_sum" # names of variables to test
)
method1 <- "anova" # one of "anova" or "kruskal.test"
method2 <- "t.test" # one of "wilcox.test" or "t.test"
my_comparisons <- list(c("EUR", "AFR"),
                       c("EUR", "EAS"),
                       c("EUR", "AMR"),
                       c("EUR", "SAS"),
                       c("EAS", "AMR"),
                       c("AFR", "AMR"),
                       c("SAS", "AFR"),
                       c("SAS", "EAS")
                       ) # comparisons for post-hoc tests
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

#ggsave("1000GP_populations_12k_ANOVA.png", width = 30, height = 13, units = "cm",dpi = "retina")


dat$spp <- factor(dat$spp, levels=c("EUR", "AFR", "AMR","EAS","SAS"))


# Perform the test
gugu<-compare_means(effect_sum ~ spp,  data = dat, p.adjust.method = "bonferroni",
              ref.group = "EUR", method = "t.test")
#write.csv(gugu,"t-test_1000g_racial.csv",row.names = TRUE)

bara<-ggboxplot(dat,
                   x = "spp", y = "effect_sum",
                   color = "spp",
                   legend = "none",
                   palette = "npg",
                   add = "jitter"
                )+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(dat$effect_sum), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 300)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "EUR")+
xlab('Population') +
  ylab('PGS')

bara<- bara + theme_cowplot(font_size = 7,line_size = 1)+ theme(legend.position = "none")
bara


ggsave("~/Desktop/Figure3_mod.pdf", bara, width=17.5, height=10, units="cm", dpi=300)

#ggsave("1000GP_populations_12k_ANOVA_t-test_EU_base.png", width = 30, height = 13, units = "cm",dpi = "retina")


# Compute the analysis of variance
res.aov <- aov(effect_sum ~ spp, data = GP_all_2)
# Summary of the analysis
summary(res.aov)

lal<-TukeyHSD(res.aov)
lal
plot(lal)
bogo<-as.data.frame(lal$spp)
#write.csv(bogo,"tukey_12k_1000gp.csv",row.names = TRUE)

#plotting more
par(mar=c(3,8,3,3))
plot(TukeyHSD(res.aov,conf.level = 0.95),las=2,tcl = -.5)


bublic<-pairwise_comparisons(
  data = GP_all_2,
  x = spp ,
  y = effect_sum,
  type = "parametric",
  var.equal = TRUE,
  paired = FALSE,
  p.adjust.method = "bonferroni"
)

pairwise_comparisons(
  data = GP_all_2,
  x = spp ,
  y = effect_sum,
  type = "parametric",
  var.equal = FALSE,
  paired = FALSE,
  p.adjust.method = "bonferroni"
)



