library(phyloseq)
library(microbiomeMarker)
library(ggplot2)


# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

ps <- tax_glom(psm_perio, taxrank = "Genus")
ps <- subset_samples(ps, Site == "Subgingival")

mm_lefse <- run_lefse(
  ps,
  wilcoxon_cutoff = 0.05,
  group = "Sex",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2)

CladoPerio <- plot_cladogram(mm_lefse, 
               color = c("#F88B9D", "#8ECEFD")) + 
                  ggtitle("Periodontitis")+
  theme(legend.position = "bottom")

CladoPerio



# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

ps <- tax_glom(psm_healthy, taxrank = "Genus")
ps <- subset_samples(ps, Site == "Subgingival")

mm_lefse <- run_lefse(
  ps,
  wilcoxon_cutoff = 0.05,
  group = "Sex",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2)

CladoHealthy <- plot_cladogram(mm_lefse, 
               color = c("#F88B9D", "#8ECEFD")) + 
                  ggtitle("Healthy")+
  theme(legend.position = "bottom")
CladoHealthy

library(grid)
pdf("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Cladograms.pdf",
    height = 9, width = 20, compress = FALSE)
 grid.draw(cbind(ggplotGrob(CladoHealthy), ggplotGrob(CladoPerio)))
dev.off()



