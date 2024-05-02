library(phyloseq)
library(ggplot2)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PlyloSeq Object for analysys.RData")


TopNOTUs <- names(sort(taxa_sums(ps), TRUE)[1:15])
ent10   <- prune_taxa(TopNOTUs, ps)

plot_bar(ent10, "Sex", fill="Periodontitis", facet_grid = Periodontitis ~ Genus)+
  geom_bar(aes(color=Sex, fill=Sex), stat="identity", position="stack")
