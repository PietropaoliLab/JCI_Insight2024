library(phyloseq)
library(here)
library(dplyr)
library(microbiomeMarker)
library(ggplot2)
library(ggpubr)
load(file = here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Syntetic_16S_dataset.RData"))


# my results df
results <- data.frame()
disease <- unique(ps@sam_data$Perio)
smoke <-  unique(ps@sam_data$Smoking)

for (d in disease) {
  for (s in smoke) {
    
    # Define my subseto of interest
    # d = "Periodontitis"
    # s = "No"
    tmp <- subset_samples(ps, Perio == d & Smoking == s)
    
    # Run t-test for DA
    foo <- run_simple_stat(tmp,
                           taxa_rank = "Genus",
                           transform = "log10p",
                           norm = "TSS", #  total sum scaling, also referred to as "relative abundance", the abundances were normalized by dividing the corresponding sample library size.
                           method = "t.test",
                           # p_adjust = "fdr",
                           group = "Sex")
    
    res <- data.frame(marker_table(foo))
    
    if (nrow(res) <1 ) {
      
      res <- data.frame("feature" = NA, 
                        "enrich_group" = NA,
                        "ef_diff_mean" = NA,
                        "pvalue"= NA,
                        "padj" = NA,
                        "Perio" = NA,
                        "Smoking" = NA) 
      
    }
    res$Perio <- d
    res$Smoking <- s
    
    results <- rbind(results, res)
    
  }
}

rm(d, foo, res, s, tmp)
# Apply FDR
results$padj <- p.adjust(results$pvalue, method = "fdr")


####################################
#
##  Microgenderome Index - MGI
#
####################################

MGI_Alpha <- data.frame()

for (d in disease) {
  for (s in smoke) {
    
# Define my subseto of interest
# d = "Healthy"
# s = "Yes"
tmp <- subset_samples(ps, Perio == d & Smoking == s)
tmp <- tax_glom(tmp, taxrank = "Genus")
tmp <- microbiome::transform(tmp, transform = "log10p")
tmp <- microbiome::transform(tmp, transform = "compositional")


# Subset Genera of interest
foo <- subset(results, Perio == d & Smoking == s)
MaleGenus <- as.character(na.omit(ifelse(foo$ef_diff_mean < 0, foo$feature, NA)))
FemaleGenus <- as.character(na.omit(ifelse(foo$ef_diff_mean > 0, foo$feature, NA)))

# Prepare data for MGI
df <- subset_taxa(tmp, Genus %in% foo$feature)
df <- otu_table(df)
df <- as.data.frame(t(df))
colnames(df) <- gsub(".*Genus_", "", colnames(df))

# Sums relative abundance of enriched genus in male or in female
df$MaleEnriched <- rowSums(df[MaleGenus])
df$FemaleEnriched <- rowSums(df[FemaleGenus])

# Compute MGI
MGI <- log((1+df$FemaleEnriched)/(1+df$MaleEnriched))


# Computing Alpha diversity using Observed
AlphaDivesity <- plot_richness(subset_samples(ps, Perio == d & Smoking == s), x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
AlphaDivesity <- AlphaDivesity[["data"]]


x <- data.frame(AlphaDivesity, MGI)

MGI_Alpha <- rbind(MGI_Alpha, x)
rm(x)

  }
}

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
ggscatter(MGI_Alpha, x = "AgeYrs", y = "MGI", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          palette = MyPalette,
          shape = 21,
          size = 3,
          fill = "Sex")+
  #add.params = list(color = "blue", fill = "lightgray"))+
  labs(x="Age (years)", y="MSI", title = "")+
  scale_x_continuous(breaks = seq(30, 75, by = 10))+
  MyTheme+
  stat_cor(aes(color=Sex),method = "pearson")+
  facet_grid(Perio ~ Smoking)


ggscatter(MGI_Alpha, x = "MGI", y = "value", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          palette = MyPalette,
          shape = 21,
          size = 3,
          fill = "Sex")+
  #add.params = list(color = "blue", fill = "lightgray"))+
  labs(x="MSI", y="Observed ASV diversity", title = "")+
  # scale_x_continuous()+
  MyTheme+
  stat_cor(aes(color=Sex),method = "pearson")+
  facet_wrap(Perio ~ Smoking, scales = "free")


ggscatter(MGI_Alpha, x = "MGI", y = "value", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          palette = MyPalette,
          shape = 21,
          size = 3,
          fill = "Sex")+
  #add.params = list(color = "blue", fill = "lightgray"))+
  labs(x="MSI", y="Observed ASV diversity", title = "")+
  # scale_x_continuous()+
  MyTheme+
  stat_cor(aes(color=Sex),method = "pearson")+
  facet_wrap(Perio ~ Smoking, scales = "free")


ggscatter(MGI_Alpha, x = "MGI", y = "value",
                              add = "reg.line",                                       # Add regression line
                              conf.int = TRUE,                                        # Add confidence interval
                              fullrange = TRUE,
                              shape = 21,
                              palette = c("#F88B9D", "#8ECEFD"),
                              alpha = 0.6,
                              size = 3,
                              fill = "Sex",
                              color = "Sex",
                              add.params = list(color = "#252F63", fill = "gray90", size = .5))+
  labs(x="MGI (Log)", y="Observed ASV diversity", title = "")+
  MyTheme+
  # coord_cartesian(ylim = c(0, 360))+
  # scale_y_continuous(breaks = seq(0, 360, by=10))+
  stat_cor(method = "pearson")+
  facet_grid(Perio ~ Smoking)




ggplot(MGI_Alpha, aes(x=Sex, y=MGI, color = Sex))+             
  geom_violin(aes(fill=Sex), alpha=.3)+
  geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=3, stroke = .7)+
  geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
  labs(title = "Healthy", y = "MGI", x="")+
  #scale_y_continuous(breaks = c(20, 40, 60, 80))+
  #coord_fixed(ylim = c(10,90))+
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme+
  theme(aspect.ratio = 2/1)+
  stat_compare_means(comparisons = list(c("female", "male")), method = "wilcox.test", label = "p.signif", family="Avenir")+
  facet_grid(Perio ~ Smoking)

