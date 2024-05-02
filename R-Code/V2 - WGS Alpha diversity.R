library(microbiomeSeq)  #load the package
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library(patchwork)
library(here)
library(readr)
loadfonts(device = "pdf")

# Load data relative to diversity
alpha <- read.csv(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS diversity_analysis/AlphaDiversity.tsv"), sep="")
Shannon <- read.csv(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS diversity_analysis/AlphaDiversity_Shannon.tsv"), sep="")
beta <- read.csv(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS diversity_analysis/BetaDiversity.tsv"), sep="")

# Assigning names to merge
alpha$Run <- gsub("_metagenome", "", rownames(alpha))
Shannon$Run <- gsub("_metagenome", "", rownames(Shannon))
rownames(beta) <- gsub("_metagenome", "", rownames(beta))

# Load Phyloseq
load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS_for_revision.RData"))


# Append observed Alpha diversity to Metadata of PhyloSeq
alpha <- merge(data.frame(wgs@sam_data), alpha)
Shannon <- merge(data.frame(wgs@sam_data), Shannon)


# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 2.5/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))



alpha$Sex <- factor(alpha$Sex, levels = c("female", "male"), labels = c("F", "M"))
Shannon$Sex <- factor(Shannon$Sex, levels = c("female", "male"), labels = c("F", "M"))


ggplot(alpha, aes(x = Sex, y = observed)) + 
geom_violin(aes(fill=Sex), alpha=.3)+
  geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=2.4)+
  geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
  labs(title = "Observed diversity", y = "Observed diversity")+
  scale_y_continuous(breaks = c(20, 40, 60, 80))+
  coord_fixed(ylim = c(10,250))+
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  facet_wrap(Smoke ~ Disease)+
  MyTheme+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family="Avenir")



ggplot(Shannon, aes(x = Sex, y = diversity_shannon)) + 
  geom_violin(aes(fill=Sex), alpha=.3)+
  geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=2.4)+
  geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
  labs(title = "Shannon diversity", y = "Observed diversity")+
  scale_y_continuous(breaks = c(2, 4, 6))+
  coord_fixed(ylim = c(0,6))+
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  facet_wrap(Smoke ~ Disease)+
  MyTheme+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family="Avenir")



