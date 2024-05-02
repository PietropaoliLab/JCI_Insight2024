library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ANCOMBC)
library(dplyr)
library(DESeq2)
library(microbiomeMarker)

set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

MySite <- unique(psm_perio@sam_data$Site)

# LEfSe Periodontitis
lefse <- data.frame()
for (i in MySite) {
foo <- run_lefse(subset_samples(pseq_genus, Site == i), 
                         group = "Sex",
                         taxa_rank = "Genus",
                         norm = "TSS", # normalized the data using total sum scaling
                         kw_cutoff = 1, # Then it performed a welch.test-Wallis 
                         wilcoxon_cutoff = 1,
                         lda_cutoff = 2) # threshold score of 2.0 (default)
foo <- data.frame(foo@marker_table)
foo$Genus <- foo$feature
foo$Model <- "LEfSe"
foo$Site <- i

lefse <- rbind(lefse, foo)
rm(foo, i)
}


# ANCOM-BC Periodontitis
ancom <- data.frame()
for (i in MySite) {
foo <- ancombc(phyloseq = subset_samples(pseq_genus, Site == i),
                 formula = "Sex", 
                 p_adj_method = "fdr", 
                 zero_cut = 0.90, # by default prevalence filter of 10% is applied
                 lib_cut = 0, 
                 group = "Sex", 
                 struc_zero = TRUE, 
                 neg_lb = TRUE, 
                 tol = 1e-5, 
                 max_iter = 100, 
                 conserve = TRUE, 
                 alpha = 0.05, 
                 global = TRUE)

foo <- data.frame(foo$res)
colnames(foo) <- c("beta", "se", "W", "p_val", "padj", "diff_abn")
foo$Genus <- rownames(foo)
foo$Genus <- gsub(".*Genus_", "", foo$Genus)
foo$Model = "ANCOM-BC"
foo$Site = i
ancom <- rbind(ancom, foo)
rm(foo, i)
}

# DESeq2 Periodontitis
Deseq <- data.frame()
for (i in MySite) {
foo <- run_deseq2(subset_samples(pseq_genus, Site == i), 
                  group = "Sex",
                  taxa_rank = "Genus",
                  p_adjust = "fdr",
                  norm = "RLE", # relative log expression
                  sfType = "poscounts",
                  fitType='local',
                  contrast = c("male", "female"))
  
foo <- data.frame(foo@marker_table)

names(foo) <- c("Genus","enrich_group", "ef_logFC","pvalue","padj") 
foo$Model = "DESeq2"
foo$Site = i
Deseq <- rbind(Deseq, foo)
rm(foo, i)
}



# welch.test
welch <- data.frame()
for (i in MySite) {
  foo <- run_simple_stat(subset_samples(pseq_genus, Site == i),
                       taxa_rank = "Genus",
                       transform = "log10p",
                       norm = "TSS",
                       method = "welch.test",
                       # p_adjust = "fdr",
                       group = "Sex")

foo <- data.frame(marker_table(foo))
names(foo) <- c("Genus","enrich_group", "ef_diff_mean","pvalue","padj") 

foo$Model = "Welch's test"
foo$Site = i
foo <- foo[,c("Genus", "padj", "Model", "Site")]
welch <- rbind(welch, foo)
rm(foo, i)
}

# FDR correction
welch$padj <- p.adjust(welch$padj, method = "fdr")

# ALDEX2
ALDEx2 <- data.frame()
for (i in MySite) {
  foo <- run_aldex(subset_samples(pseq_genus, Site == i), 
                   taxa_rank = "Genus",
                   transform = "identity", # Raw reads is recommended from the ALDEx2 paper
                   norm = "RLE",#relative log expression,
                   p_adjust = "fdr",
                   group = "Sex")
  
  foo <- data.frame(marker_table(foo))
  names(foo) <- c("Genus","enrich_group", "ef_aldex","pvalue","padj") 
  
  foo$Model = "ALDEx2"
  foo$Site = i
  foo <- foo[,c("Genus", "padj", "Model", "Site")]
  
  ALDEx2 <- rbind(ALDEx2, foo)
  rm(foo, i)
}




# Cobine periodontitis results
res <- rbind(Deseq[,c("Genus", "padj", "Model", "Site")], ancom[,c("Genus", "padj", "Model",  "Site")])
res <- rbind(res, lefse[,c("Genus", "padj", "Model", "Site")])
res <- rbind(res, welch[,c("Genus", "padj", "Model", "Site")])
res <- rbind(res, ALDEx2[,c("Genus", "padj", "Model", "Site")])


res$Group = "Periodontitis"
table(res$Model)
table(res$Genus, res$Model)
table(res$Site, res$Model)


Perio <- res
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, i, welch, MySite, psm_perio)

################################
###       HEALTHY
################################
set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")

MySite <- unique(psm_healthy@sam_data$Site)
# LEfSe Healthy
lefse <- data.frame()
for (i in MySite) {
  foo <- run_lefse(subset_samples(pseq_genus, Site == i), 
                   group = "Sex",
                   taxa_rank = "Genus",
                   norm = "TSS", # normalized the data using total sum scaling
                   kw_cutoff = 1, # Then it performed a welch.test-Wallis 
                   wilcoxon_cutoff = 1,
                   lda_cutoff = 2) # threshold score of 2.0 (default)
  foo <- data.frame(foo@marker_table)
  foo$Genus <- foo$feature
  foo$Model <- "LEfSe"
  foo$Site <- i
  
  lefse <- rbind(lefse, foo)
  rm(foo, i)
  
}

# ANCOM-BC Healthy
ancom <- data.frame()
for (i in MySite) {
  foo <- ancombc(phyloseq = subset_samples(pseq_genus, Site == i),
                 formula = "Sex", 
                 p_adj_method = "fdr", 
                 zero_cut = 0.90, # by default prevalence filter of 10% is applied
                 lib_cut = 0, 
                 group = "Sex", 
                 struc_zero = TRUE, 
                 neg_lb = TRUE, 
                 tol = 1e-5, 
                 max_iter = 100, 
                 conserve = TRUE, 
                 alpha = 0.05, 
                 global = TRUE)
  
  foo <- data.frame(foo$res)
  colnames(foo) <- c("beta", "se", "W", "p_val", "padj", "diff_abn")
  foo$Genus <- rownames(foo)
  foo$Genus <- gsub(".*Genus_", "", foo$Genus)
  foo$Model = "ANCOM-BC"
  foo$Site = i
  ancom <- rbind(ancom, foo)
  rm(foo, i)
  
}


Deseq <- data.frame()
# DESeq2 Healthy
for (i in MySite) {
  
  foo <- run_deseq2(subset_samples(pseq_genus, Site == i), 
                    group = "Sex",
                    taxa_rank = "Genus",
                    p_adjust = "fdr",
                    norm = "RLE", # relative log expression
                    sfType = "poscounts",
                    fitType='local',
                    contrast = c("male", "female"))
  
  foo <- data.frame(foo@marker_table)
  
  names(foo) <- c("Genus","enrich_group", "ef_logFC","pvalue","padj") 
  foo$Model = "DESeq2"
  foo$Site = i
  
  
  Deseq <- rbind(Deseq, foo)
  rm(foo, i)
  
}

# welch.test Healthy
welch <- data.frame()
for (i in MySite) {
  
  
  foo <- run_simple_stat(subset_samples(pseq_genus, Site == i),
                         taxa_rank = "Genus",
                         transform = "log10p",
                         norm = "TSS",
                         method = "welch.test",
                         # p_adjust = "fdr",
                         group = "Sex")
  
  foo <- data.frame(marker_table(foo))
  names(foo) <- c("Genus","enrich_group", "ef_diff_mean","pvalue","padj") 
  
  foo$Model = "Welch's test"
  foo$Site = i
  foo <- foo[,c("Genus", "padj", "Model", "Site")]
  
  welch <- rbind(welch, foo)
  rm(foo, i)
}

# FDR correction
welch$padj <- p.adjust(welch$padj, method = "fdr")


# ALDEX2 Healthy
ALDEx2 <- data.frame()
for (i in MySite) {
  foo <- run_aldex(subset_samples(pseq_genus, Site == i), 
                   taxa_rank = "Genus",
                   transform = "identity", # Raw reads is recommended from the ALDEx2 paper
                   norm = "RLE",#relative log expression,
                   p_adjust = "fdr",
                   group = "Sex")
  
  
  foo <- data.frame(marker_table(foo))
  names(foo) <- c("Genus","enrich_group", "ef_aldex","pvalue","padj") 
  
  foo$Model = "ALDEx2"
  foo$Site = i
  foo <- foo[,c("Genus", "padj", "Model", "Site")]
  
  ALDEx2 <- rbind(ALDEx2, foo)
  rm(foo, i)
}




# Cobine Healthy results
res <- rbind(Deseq[,c("Genus", "padj", "Model", "Site")], ancom[,c("Genus", "padj", "Model",  "Site")])
res <- rbind(res, lefse[,c("Genus", "padj", "Model", "Site")])
res <- rbind(res, welch[,c("Genus", "padj", "Model", "Site")])
res <- rbind(res, ALDEx2[,c("Genus", "padj", "Model", "Site")])


res$Group = "Healthy"
table(res$Model)
table(res$Genus, res$Model)
table(res$Site, res$Model)

Healthy <- res


results <- rbind(Perio, Healthy)

rownames(results) <- seq(1:nrow(results))
results$SigGenera <- ifelse(results$padj <0.05, 1, NA)


rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_healthy, res)




library(ggplot2)
results$Model <- factor(results$Model)
Periodontitis <- subset(results, Group == "Periodontitis" & SigGenera == 1)
Healthy <- subset(results, Group == "Healthy"  & SigGenera == 1)


MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8.5),
        #panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 2/.6,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))

p1 <- ggplot(Periodontitis) +
  geom_point(aes(x = Model, y = Genus, fill=SigGenera), fill = "deeppink", size = 3, shape=21)+
  facet_wrap(~Site, scales="free_y")+
  labs(x="", y="", title = "Periodontitis")+
  MyTheme


p2 <- ggplot(Healthy) +
  geom_point(aes(x = Model, y = Genus, fill=SigGenera), fill = "deeppink", size = 3, shape=21)+
  scale_x_discrete(drop=FALSE)+ # keep ALDEX2 on
  facet_wrap(~Site, scales="free_y", nrow = 1)+
  labs(x="", y="", title = "Healthy")+
  MyTheme

library(gridExtra)
grid.arrange(p1, p2, nrow = 2)

p1
p2



ggplot(Periodontitis) +
  geom_point(aes(x = Model, y = Genus, fill = padj), size = 4, shape=21)+
  scale_fill_gradient2(low = "mediumseagreen",
                       mid = "yellow",
                       high = "#6495ED",
                       midpoint = 0.025)+
  facet_wrap(~Site, scales="free_y")+
  labs(x="", y="", title = "Periodontitis")+
  MyTheme


ggplot(Healthy) +
  geom_point(aes(x = Model, y = Genus, fill = padj), size = 4, shape=21)+
  scale_fill_gradient2(low = "mediumseagreen",
                       mid = "yellow",
                       high = "#6495ED",
                       midpoint = 0.025)+
  facet_wrap(~Site, scales="free_y")+
  labs(x="", y="", title = "Healthy")+
  MyTheme


######################################à
#
## Defining Genus for Genderome Index
#
######################################à

library(plyr)
# Genus of Interest
goi <- ddply(results, c("Genus", "Site", "Group"), summarise,
             ModelsSigns = sum(SigGenera, na.rm = TRUE))

# subset only when 3 or more method discover genus
goi <- subset(goi, ModelsSigns >= 3)

# subsetting only subgingival
goi <- subset(goi, Site == "Subgingival" & Group != "Healthy")

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Subset only subgingival site
pseq_genus <- subset_samples(pseq_genus, Site == "Subgingival")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")
taxa_names(pseq_genus)

# Subsetting genus of interest
pseq_genus <- subset_taxa(pseq_genus, Genus %in% goi$Genus)
taxa_names(pseq_genus)

# Using identity
pseq_genus <- microbiome::transform(pseq_genus, transform = "log10p")

goi_plot <- psmelt(pseq_genus)




library(ggpubr)
ggboxplot(goi_plot, x = "Sex", y = "Abundance",
          add = "mean",
          facet.by = "Genus",
          scales = "free")+
  stat_compare_means()+
  stat_summary(fun.data = function(x) data.frame(y=5, label = paste(round(mean(x), 3))), geom="text") 
  









