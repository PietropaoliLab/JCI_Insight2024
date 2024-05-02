library(phyloseq)
library(ggpubr)
library(ggplot2)
library(microbiome)
library(plyr)
library(reshape2)
library(tidyr)


# Load my 16S data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")


res <- data.frame()
#################
# Periodontitis
################
bioproject <- unique(psm_perio@sam_data$BioProject)

for (b in bioproject) {
  #  b="PRJEB6047"
  tmp <- psm_perio
  tmp <- tax_glom(tmp, taxrank = "Phylum")
  tmp <- subset_samples(tmp, BioProject == b)
  site <- unique(tmp@sam_data$Site)
  
  # Clycle for Perio
  for (s in site) {
     # s="Plaque"
    cat(paste0("Periodontitis - ", b, ": ", s, "\n"))
    foo <- phyloseq::subset_samples(tmp, Site == s)
    # Relative abundances
    foo <- microbiome::transform(foo, "compositional")
    foo <- psmelt(foo)
    foo$Sex <- factor(foo$Sex, levels = c("female", "male"), labels = c("F", "M"))
    foo <- subset(foo, Abundance >0.01)
    
    # compute FDR with wilcoxon
    fdr <- compare_means(Abundance ~ Sex, group.by = "Phylum", method = "wilcox.test", p.adjust.method = "fdr", data = foo)
    
    # compute average between gender
    x <- ddply(foo, c("Sex", "Phylum"), summarize,
               avg = mean(Abundance, na.rm = TRUE),
               N = sum(!is.na(Abundance)))
    x <- dcast(x, Phylum ~ Sex, value.var = "avg")
    
    fdr <- merge(fdr, x)
    fdr$BioProject <- b
    fdr$Site <- s
    fdr$Periodontitis <- "Periodontitis"
    
    res <- rbind(res, fdr)
  }
}

###################
# Healthy
##################
bioproject <- unique(psm_healthy@sam_data$BioProject)

for (b in bioproject) {
  #  b="PRJEB6047"
  tmp <- psm_healthy
  tmp <- tax_glom(tmp, taxrank = "Phylum")
  tmp <- subset_samples(tmp, BioProject == b)
  site <- unique(tmp@sam_data$Site)
  
  # Clycle for Perio
  for (s in site) {
    # s="Plaque"
    cat(paste0("Healthy - ", b, ": ", s, "\n"))
    foo <- phyloseq::subset_samples(tmp, Site == s)
    # Relative abundances
    foo <- microbiome::transform(foo, "compositional")
    foo <- psmelt(foo)
    foo$Sex <- factor(foo$Sex, levels = c("female", "male"), labels = c("F", "M"))
    foo <- subset(foo, Abundance >0.01)
    
    # compute FDR with wilcoxon
    fdr <- compare_means(Abundance ~ Sex, group.by = "Phylum", method = "wilcox.test", p.adjust.method = "fdr", data = foo)
    
    # compute average between gender
    x <- ddply(foo, c("Sex", "Phylum"), summarize,
               avg = mean(Abundance, na.rm = TRUE),
               N = sum(!is.na(Abundance)))
    x <- dcast(x, Phylum ~ Sex, value.var = "avg")
    
    fdr <- merge(fdr, x)
    fdr$BioProject <- b
    fdr$Site <- s
    fdr$Periodontitis <- "Healthy"
    
    res <- rbind(res, fdr)
  }
}



res$diff <- res$F - res$M
res$gFC <- log(res$F) - log(res$M)

df <- melt(res, id.vars = c("Phylum", "Site", "Periodontitis", "BioProject"), measure.vars = c("F", "M", "p.adj"))
names(df) <- c("Phylum", "Site", "Periodontitis", "Sex", "Abundance")




MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 2/0.8,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(fill = NA, colour = NA))


p1 <- ggplot(Healthy, aes(x=Sex, y=Phylum, fill = Abundance))+
  geom_tile(color = "gray15")+
  facet_grid(~ Site)+
  scale_fill_distiller(palette = "BuGn", direction = 1)+
  MyTheme +
  labs(x="", y="", title = "Healthy")+
  theme(plot.margin = unit(c(0,0,-2,-2), "mm"),
        panel.border = element_rect(colour = NA),
        panel.spacing = unit(c(-2,-2), "mm"))


p2 <- ggplot(Perio, aes(x=Sex, y=Phylum, fill = Abundance))+
  geom_tile(color = "gray15")+
  facet_grid(~ Site)+
  scale_fill_distiller(palette = "BuGn", direction = 1)+
  MyTheme +
  labs(x="", y="", title = "Periodontitis")+
  theme(plot.margin = unit(c(0,0,-2,-2), "mm"),
        panel.border = element_rect(colour = NA),
        panel.spacing = unit(c(-2,-2), "mm"))

ggarrange(p1, p2)

ggsave(filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Heatmaps.eps",
       plot = last_plot(), width = 6, height = 5.5)
