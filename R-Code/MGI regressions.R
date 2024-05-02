library(phyloseq)
library(dplyr)
library(microbiomeMarker)
library(ggpubr)
library(ggplot2)
library(here)

load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData"))

# my results df
results <- data.frame()

# Define my subseto of interest
tmp <- subset_samples(psm_perio, Site == "Subgingival" & Smoking == "No")

# Run t-test for DA
foo <- run_simple_stat(tmp,
                       taxa_rank = "Genus",
                       transform = "log10p",
                       norm = "TSS", #  total sum scaling, also referred to as "relative abundance", the abundances were normalized by dividing the corresponding sample library size.
                       method = "t.test",
                       # p_adjust = "fdr",
                       group = "Sex")

results <- data.frame(marker_table(foo))
results$Model = "t-test"
results$Smoking = "Non-smoking"

# Apply FDR
results$padj <- p.adjust(results$pvalue, method = "fdr")

# Plotting results
ggplot(results, aes(x=reorder(feature, ef_diff_mean), y = ef_diff_mean, fill = enrich_group))+
  geom_col()+
  facet_wrap(~Smoking, scales = "free")+
  coord_flip()

####################################
#
##  Microgenderome Index - MGI
#
####################################

tmp <- tax_glom(tmp, taxrank = "Genus")
tmp <- microbiome::transform(tmp, transform = "log10p")
tmp <- microbiome::transform(tmp, transform = "compositional")


# Subset Genera of interest
MaleGenus <- na.omit(ifelse(results$ef_diff_mean < 0, results$feature, NA))
FemaleGenus <- na.omit(ifelse(results$ef_diff_mean > 0, results$feature, NA))

# Prepare data for MGI
df <- subset_taxa(tmp, Genus %in% results$feature)
df <- data.frame(otu_table(df))
colnames(df) <- gsub(".*Genus_", "", colnames(df))

# Sums relative abundance of enriched genus in male or in female
df$MaleEnriched <- rowSums(df[MaleGenus])
df$FemaleEnriched <- rowSums(df[FemaleGenus])

# Compute MGI
df$MGI <- log((1+df$FemaleEnriched)/(1+df$MaleEnriched))

# Store MGI into tmp phylobject
tmp@sam_data$MGI <- df$MGI

# Computing Alpha diversity using Observed
alpha <- plot_richness(tmp, x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
df <- alpha[["data"]]

# Plotting results
MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))

MyPalette <- c("#F88B9D", "#8ECEFD")
ggscatter(df, x = "AgeYrs", y = "MGI", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          palette = MyPalette,
          shape = 21,
          size = 3,
          fill = "Sex")+
          #add.params = list(color = "blue", fill = "lightgray"))+
  labs(x="Age (years)", y="MGI", title = "Periodontitis")+
  scale_x_continuous(breaks = seq(30, 75, by = 10))+
  MyTheme+
  stat_cor(aes(color=Sex),method = "pearson")


# Observed asv regression in periodontitis
RegPerioVsAge <- ggscatter(df, x = "AgeYrs", y = "value", color = "Sex",
                          add = "reg.line",                                       # Add regression line
                          conf.int = TRUE,                                        # Add confidence interval
                          fullrange = TRUE,
                          palette = MyPalette,
                          shape = 21,
                          size = 3,
                          fill = "Sex")+
                          #add.params = list(color = "blue", fill = "lightgray"))+
                          labs(x="Age (years)", y="Observed ASV diversity", title = "Periodontitis")+
                          scale_x_continuous(breaks = seq(30, 75, by = 10))+
                          MyTheme+
                          stat_cor(aes(color=Sex),method = "pearson")



########################################################
#
#
#  ***************    HEALTHY     ***************       #
#
########################################################


library(phyloseq)
library(dplyr)
library(microbiomeMarker)
library(ggpubr)
library(ggplot2)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# my results df
results <- data.frame()

# Define my subseto of interest
tmp <- subset_samples(psm_healthy, Site == "Subgingival" & Smoking == "No")

# Run t-test for DA
foo <- run_simple_stat(tmp,
                       taxa_rank = "Genus",
                       transform = "log10p",
                       norm = "TSS", #  total sum scaling, also referred to as "relative abundance", the abundances were normalized by dividing the corresponding sample library size.
                       method = "t.test",
                       # p_adjust = "fdr",
                       group = "Sex")

results <- data.frame(marker_table(foo))
results$Model = "t-test"
results$Smoking = "Non-smoking"

# Apply FDR
results$padj <- p.adjust(results$pvalue, method = "fdr")

# Plotting results
ggplot(results, aes(x=reorder(feature, ef_diff_mean), y = ef_diff_mean, fill = enrich_group))+
  geom_col()+
  facet_wrap(~Smoking, scales = "free")+
  coord_flip()

####################################
#
##  Microgenderome Index - MGI
#
####################################

tmp <- tax_glom(tmp, taxrank = "Genus")
tmp <- microbiome::transform(tmp, transform = "log10p")
tmp <- microbiome::transform(tmp, transform = "compositional")


# Subset Genera of interest
MaleGenus <- na.omit(ifelse(results$ef_diff_mean < 0, results$feature, NA))
FemaleGenus <- na.omit(ifelse(results$ef_diff_mean > 0, results$feature, NA))

# Prepare data for MGI
df <- subset_taxa(tmp, Genus %in% results$feature)
df <- data.frame(otu_table(df))
colnames(df) <- gsub(".*Genus_", "", colnames(df))

# Sums relative abundance of enriched genus in male or in female
df$MaleEnriched <- rowSums(df[MaleGenus])
df$FemaleEnriched <- rowSums(df[FemaleGenus])

# Compute MGI
df$MGI <- log((1+df$FemaleEnriched)/(1+df$MaleEnriched))

# Store MGI into tmp phylobject
tmp@sam_data$MGI <- df$MGI

# Computing Alpha diversity using Observed
alpha <- plot_richness(tmp, x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
df <- alpha[["data"]]

# Plotting results
MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))


ggscatter(df, x = "AgeYrs", y = "MGI", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          shape = 21,
          size = 3,
          fill = "Sex")+
  #add.params = list(color = "blue", fill = "lightgray"))+
  labs(x="Age (years)", y="MGI", title = "Healthy")+
  MyTheme+
  stat_cor(aes(color=Sex),method = "pearson")


# Shannon
RegHealthyVsAge <- ggscatter(df, x = "AgeYrs", y = "value", color = "Sex",
                          add = "reg.line",                                       # Add regression line
                          conf.int = TRUE,                                        # Add confidence interval
                          fullrange = TRUE,
                          palette = MyPalette,
                          shape = 21,
                          size = 3,
                          fill = "Sex")+
                          #add.params = list(color = "blue", fill = "lightgray"))+
                          labs(x="Age (years)", y="Observed ASV diversity", title = "Healthy")+
                          MyTheme+
                          stat_cor(aes(color=Sex),method = "pearson")




RegHealthyVsAge + RegPerioVsAge