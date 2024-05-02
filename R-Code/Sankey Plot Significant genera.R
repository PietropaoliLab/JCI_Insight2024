library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ANCOMBC)
library(dplyr)
library(plyr)
library(DESeq2)
library(microbiomeMarker)
library(ggplot2)
library(extrafont)
library(showtext)
library(patchwork)
library(viridis)

loadfonts(device = "pdf")

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
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_perio)







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

# store my results
results <- rbind(Perio, Healthy)

# changing row name and identify my genera with FDR <0.05
rownames(results) <- seq(1:nrow(results))
results$SigGenera <- ifelse(results$padj <0.05, 1, NA)

# cleaning
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_healthy, res)


# Subsetting my data
results$Model <- factor(results$Model)

# Ordering Sampling Sites
results$Site <- factor(results$Site, levels = c("Saliva", "Plaque", "Subgingival"))

# Subsetting my dataset for healthy and periodontitis 
Periodontitis <- subset(results, Group == "Periodontitis" & SigGenera == 1)
Healthy <- subset(results, Group == "Healthy"  & SigGenera == 1)

# Overall database
Overall <- rbind(Periodontitis, Healthy)
Overall$Smoking <- "Overall"


######################## SMOKERS START HERE#  ######################################
#
#                           SMOKERS

set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")
pseq_genus <- subset_samples(pseq_genus, Smoking == "Yes") # Subletting only smoker individuals
MySite <- unique(pseq_genus@sam_data$Site)

###      ---> Periodontitis <---
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
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_perio)




###     -------   HEALTHY

set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")
pseq_genus <- subset_samples(pseq_genus, Smoking == "Yes")


MySite <- unique(pseq_genus@sam_data$Site)
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

# store my results
results <- rbind(Perio, Healthy)

# changing row name and identify my genera with FDR <0.05
rownames(results) <- seq(1:nrow(results))
results$SigGenera <- ifelse(results$padj <0.05, 1, NA)

# cleaning
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_healthy, res)


# Subsetting my data
results$Model <- factor(results$Model)

# Ordering Sampling Sites
results$Site <- factor(results$Site, levels = c("Saliva", "Plaque", "Subgingival"))

# Subsetting my dataset for healthy and periodontitis 
Periodontitis <- subset(results, Group == "Periodontitis" & SigGenera == 1)
Healthy <- subset(results, Group == "Healthy"  & SigGenera == 1)

# Smokers dataframe
Smokers <- rbind(Periodontitis, Healthy)
Smokers$Smoking <- "Smokers"

######################## SMOKER END HERE  ########################################


######################## !!!!! NON SMOKERS START HERE#  ######################################
set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")
pseq_genus <- subset_samples(pseq_genus, Smoking == "No") # Subletting only smoker individuals
MySite <- unique(pseq_genus@sam_data$Site)

##--------------> Periodontitis <---


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
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_perio)





#------------------> HEALTHY

set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")
pseq_genus <- subset_samples(pseq_genus, Smoking == "No")


MySite <- unique(pseq_genus@sam_data$Site)
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

# store my results
results <- rbind(Perio, Healthy)

# changing row name and identify my genera with FDR <0.05
rownames(results) <- seq(1:nrow(results))
results$SigGenera <- ifelse(results$padj <0.05, 1, NA)

# cleaning
rm(ancom, Deseq, ALDEx2, lefse, pseq_genus, welch, MySite, psm_healthy, res)


# Subsetting my data
results$Model <- factor(results$Model)

# Ordering Sampling Sites
results$Site <- factor(results$Site, levels = c("Saliva", "Plaque", "Subgingival"))

# Subsetting my dataset for healthy and periodontitis 
Periodontitis <- subset(results, Group == "Periodontitis" & SigGenera == 1)
Healthy <- subset(results, Group == "Healthy"  & SigGenera == 1)

# Smokers dataframe
NonSmokers <- rbind(Periodontitis, Healthy)
NonSmokers$Smoking <- "Non Smokers"

######################## NON SMOKER END HERE  ########################################

df <- rbind(Smokers, NonSmokers)

# Compute number of agreement
m <- ddply(df, c("Genus", "Group","Site", "Smoking"), summarise,
           Agreement = sum(!is.na(SigGenera)))
df <- merge(df, m, all.x = TRUE)

Sankey <- subset(df, Agreement > 2)

Sankey <- Sankey[,-c(5:7)]
CombinationAnalysis <- Sankey

Sankey <- unique(Sankey)

# xlsx::write.xlsx(Sankey, file = "AQ Projects/Microbiome/Gender features of periodontal microbiome/Sankey data.xlsx")


library(xlsx)
library(easyalluvial)
library(reshape2)
library(tidyverse)
library(plyr)

Sankey <- read.xlsx(file = "AQ Projects/Microbiome/Gender features of periodontal microbiome/Sankey data.xlsx", sheetIndex = 1)
Sankey <- Sankey[-1]


col_vector <- viridis::viridis(option = "C", 28)

Sankey <- Sankey[,c(1,3,4,2)]

p <- alluvial_wide(Sankey, 
                   stratum_width = 1/8,
                   #col_vector_flow = col_vector,
                   fill_by = 'first_variable')+
  theme_void()

p



# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
Saliva <- subset(Sankey, Site == "Saliva" & Smoking == "Smokers")$Genus
Plaque <- subset(Sankey, Site == "Plaque" & Smoking == "Smokers")$Genus
Subgingival <- subset(Sankey, Site == "Subgingival" & Smoking == "Smokers")$Genus

# Chart
venn.diagram(
  x = list(Saliva, Plaque, Subgingival),
  category.names = c("Saliva" , "Plaque" , "Subgingival"),
  filename = 'AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Smokers_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)





# Generate 3 sets of 200 words
Saliva <- subset(Sankey, Site == "Saliva" & Smoking == "Non Smokers")$Genus
Plaque <- subset(Sankey, Site == "Plaque" & Smoking == "Non Smokers")$Genus
Subgingival <- subset(Sankey, Site == "Subgingival" & Smoking == "Non Smokers")$Genus

# Chart
venn.diagram(
  x = list(Saliva, Plaque, Subgingival),
  category.names = c("Saliva" , "Plaque" , "Subgingival"),
  filename = 'AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Non smokers_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)


