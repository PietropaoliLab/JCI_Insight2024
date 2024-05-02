library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM PlyloSeq Object for analysys.RData")

myRun <- c(psm_healthy@sam_data$Run, psm_perio@sam_data$Run)

ps <- subset_samples(psm_ps, Run %in% myRun)
rm(psm_healthy, psm_perio, psm_ps)

ps@sam_data$Comparison <- paste(ps@sam_data$Sex, ps@sam_data$Periodontitis, sep = " ")
ps@sam_data$Comparison <- factor(ps@sam_data$Comparison, levels = c("female Healthy", "female Periodontitis", "male Healthy", "male Periodontitis"),
                                                         labels = c("FH", "FP", "MH", "MP"))
unique(ps@sam_data$Comparison)

ps1 <- prune_taxa(taxa_sums(ps) > 0, ps)
tab <- microbiome::alpha(ps1, index = "all")

ps1.meta <- meta(ps1)

ps1.meta$Shannon <- tab$diversity_shannon
# create a list of pairwise comaprisons
Comparison <- levels(ps1.meta$Comparison)

# make a pairwise list that we want to compare.
pairs <- combn(seq_along(Comparison), 2, simplify = FALSE, FUN = function(i)Comparison[i])

print(pairs)
library(ggpubr)
MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background=element_blank())


ggviolin(ps1.meta, x = "Comparison", y = "Shannon",
               add = "boxplot", fill = "Comparison", palette = c("#F88B9D", "#FC0F3A", "#8ECEFD", "#1C92D8")) +
  stat_compare_means(comparisons = pairs, method = "wilcox.test", p.adjust.methods = "bonferroni") +
  labs(x="", y="Alpha diversity (Shannon index)")+
  MyTheme +
   theme(legend.position = "none")
