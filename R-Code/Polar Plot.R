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


# ========================================================#
#                                                         # 
#       ---- >    START PLOTTING HERE     < ----          #
#                                                         #   
# ========================================================#

# Define my theme
MyTheme <- theme_bw(12)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Avenir", color = "black"),
        aspect.ratio = 1,
        legend.position = "bottom")

# Defining dimention of circle
MyLimits <- c(-1.5, 5)

# Defining my color palette
# this is a color-blindness friendly palette Ref. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
#MyPalette <-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


MyPalette <- viridis(option = "E", 5)

# Setup datasets for cycle
Site <- c("Saliva", "Plaque", "Subgingival")
MyPlots <- list()

# ===============================
# Generate plot for Periodontitis
# ===============================
for (s in Site) {
  

foo <- subset(Periodontitis, Site == s)

 # Compute number of agreement
 m <- ddply(foo, c("Genus", "Site"), summarize,
            Agreement = sum(!is.na(SigGenera)))
 foo <- merge(foo, m)
 

 
 # Get the name and the y position of each label
 label_data <- foo
 label_data <- subset(label_data[,c(1,7)], !duplicated(label_data$Genus))
 
 # Sorting data
 label_data <-label_data[order(-label_data$Agreement),]
 foo <-foo[order(-foo$Agreement),]
 
 # Creating custom id
 label_data$id <- seq(1, nrow(label_data))
 
 # get numer of bars
 number_of_bar <- nrow(label_data)
 
 # 
 angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
 label_data$hjust <- ifelse(angle < -90, 1, 0)
 label_data$angle <- ifelse(angle < -90, angle+180, angle)
 
 #foo <- merge(foo, label_data)
 
p <- ggplot(foo, aes(x=reorder(Genus, -Agreement), y=SigGenera))+
          geom_bar(stat = "identity", aes(fill = Model), width = .9)+
          scale_y_continuous(limits = MyLimits)+                # Circle dimension
          scale_fill_manual(values = MyPalette, drop = FALSE)+   # My color palette
          ggtitle("Periodontitis", subtitle = s)+                                           # My title
          coord_polar()+
          # Add labels on top of each bar
          geom_text(data = label_data, aes(x=reorder(Genus, -Agreement), y=Agreement+0.02, label=Genus, hjust=hjust), 
                    family = "Avenir", size=2.47, # 7Pt * Magic number 0.352777778
                    angle = label_data$angle, 
                    inherit.aes = FALSE) +
          MyTheme +
          theme_void()

MyPlots[[paste0("Periodontitis-", s)]] <- p

rm(foo, label_data, m, p)
  
}
rm(s)
# Store my plots
PerioPlots <- MyPlots


# reset MyPlots
MyPlots <- list()

# ===============================
# Generate plot for Healthy
# ===============================

for (s in Site) {
  
  foo <- subset(Healthy, Site == s)
  
  # Compute number of agreement
  m <- ddply(foo, c("Genus", "Site"), summarize,
             Agreement = sum(!is.na(SigGenera)))
  foo <- merge(foo, m)
  
  
  
  # Get the name and the y position of each label
  label_data <- foo
  label_data <- subset(label_data[,c(1,7)], !duplicated(label_data$Genus))
  
  # Sorting data
  label_data <-label_data[order(-label_data$Agreement),]
  foo <-foo[order(-foo$Agreement),]
  
  # Creating custom id
  label_data$id <- seq(1, nrow(label_data))
  
  # get numer of bars
  number_of_bar <- nrow(label_data)
  
  # 
  angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #foo <- merge(foo, label_data)
  
  p <- ggplot(foo, aes(x=reorder(Genus, -Agreement), y=SigGenera))+
    geom_bar(stat = "identity", aes(fill = Model), width = .9)+
    scale_y_continuous(limits = MyLimits)+                # Circle dimension
    scale_fill_manual(values = MyPalette, drop = FALSE)+   # My color palette
    ggtitle("Healthy", subtitle = s)+                                           # My title
    coord_polar()+
    # Add labels on top of each bar
    geom_text(data = label_data, aes(x=reorder(Genus, -Agreement), y=Agreement+0.02, label=Genus, hjust=hjust), 
              family = "Avenir", size=2.47, # 7Pt * Magic number 0.352777778
              angle = label_data$angle, 
              inherit.aes = FALSE) +
    MyTheme +
    theme_void()
  
  MyPlots[[paste0("Healthy-", s)]] <- p
  
  rm(foo, label_data, m, p)
  
}
rm(s)
# Store my plots
HealthyPlots <- MyPlots


# Periodontitis Plots
PerioPlots <- (PerioPlots[[1]] | PerioPlots[[2]] | PerioPlots[[3]]) 

# Healthy Plots
HealthyPlots <- (HealthyPlots[[1]] | HealthyPlots[[2]] | HealthyPlots[[3]]) 

# Combine Healthy and Periodontitis
DA_Overall <- (HealthyPlots / PerioPlots) +  plot_layout(guides = 'collect') & theme(legend.position = 'right')

DA_Overall

ggsave(plot = DA_Overall, 
       device = cairo_pdf,
       width = 11, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/DA overall.pdf")

# ================================================================================================================================


set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")
pseq_genus <- subset_samples(pseq_genus, Smoking == "Yes") # Subletting only smoker individuals

MySite <- unique(pseq_genus@sam_data$Site)

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


res$Group = "Periodontitis - Smokers"
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
                         #norm = "TSS",
                         method = "welch.test",
                         #p_adjust = "fdr",
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


# ========================================================#
#                                                         # 
#       ---- >    START PLOTTING HERE     < ----          #
#                                                         #   
# ========================================================#

# Define my theme
MyTheme <- theme_bw(12)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Avenir", color = "black"),
        aspect.ratio = 1,
        legend.position = "bottom")

# Defining dimention of circle
MyLimits <- c(-1.5, 5.5)

# Defining my color palette
# this is a color-blindness friendly palette Ref. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
#MyPalette <-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


MyPalette <- viridis(option = "E", 5)

# Setup datasets for cycle
Site <- c("Saliva", "Plaque", "Subgingival")
MyPlots <- list()

# ===============================
# Generate plot for Periodontitis
# ===============================
for (s in Site) {
  
  
  foo <- subset(Periodontitis, Site == s)
  
  # Compute number of agreement
  m <- ddply(foo, c("Genus", "Site"), summarize,
             Agreement = sum(!is.na(SigGenera)))
  foo <- merge(foo, m)
  
  
  
  # Get the name and the y position of each label
  label_data <- foo
  label_data <- subset(label_data[,c(1,7)], !duplicated(label_data$Genus))
  
  # Sorting data
  label_data <-label_data[order(-label_data$Agreement),]
  foo <-foo[order(-foo$Agreement),]
  
  # Creating custom id
  label_data$id <- seq(1, nrow(label_data))
  
  # get numer of bars
  number_of_bar <- nrow(label_data)
  
  # 
  angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #foo <- merge(foo, label_data)
  
  p <- ggplot(foo, aes(x=reorder(Genus, -Agreement), y=SigGenera))+
    geom_bar(stat = "identity", aes(fill = Model), width = .9)+
    scale_y_continuous(limits = MyLimits)+                # Circle dimension
    scale_fill_manual(values = MyPalette, drop = FALSE)+   # My color palette
    ggtitle(s)+                                           # My title
    coord_polar()+
    # Add labels on top of each bar
    geom_text(data = label_data, aes(x=reorder(Genus, -Agreement), y=Agreement+0.02, label=Genus, hjust=hjust), 
              family = "Avenir", size=2.47, # 7Pt * Magic number 0.352777778
              angle = label_data$angle, 
              inherit.aes = FALSE) +
    MyTheme +
    theme_void()
  
  MyPlots[[paste0("Periodontitis-", s)]] <- p
  
  rm(foo, label_data, m, p)
  
}
rm(s)
# Store my plots
PerioPlots <- MyPlots


# reset MyPlots
MyPlots <- list()

# ===============================
# Generate plot for Healthy
# ===============================

for (s in Site) {
  
  foo <- subset(Healthy, Site == s)
  
  # Compute number of agreement
  m <- ddply(foo, c("Genus", "Site"), summarize,
             Agreement = sum(!is.na(SigGenera)))
  foo <- merge(foo, m)
  
  
  
  # Get the name and the y position of each label
  label_data <- foo
  label_data <- subset(label_data[,c(1,7)], !duplicated(label_data$Genus))
  
  # Sorting data
  label_data <-label_data[order(-label_data$Agreement),]
  foo <-foo[order(-foo$Agreement),]
  
  # Creating custom id
  label_data$id <- seq(1, nrow(label_data))
  
  # get numer of bars
  number_of_bar <- nrow(label_data)
  
  # 
  angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #foo <- merge(foo, label_data)
  
  p <- ggplot(foo, aes(x=reorder(Genus, -Agreement), y=SigGenera))+
    geom_bar(stat = "identity", aes(fill = Model), width = .9)+
    scale_y_continuous(limits = MyLimits)+                # Circle dimension
    scale_fill_manual(values = MyPalette, drop = FALSE)+   # My color palette
    ggtitle(s)+                                           # My title
    coord_polar()+
    # Add labels on top of each bar
    geom_text(data = label_data, aes(x=reorder(Genus, -Agreement), y=Agreement+0.02, label=Genus, hjust=hjust), 
              family = "Avenir", size=2.47, # 7Pt * Magic number 0.352777778
              angle = label_data$angle, 
              inherit.aes = FALSE) +
    MyTheme +
    theme_void()
  
  MyPlots[[paste0("Healthy-", s)]] <- p
  
  rm(foo, label_data, m, p)
  
}
rm(s)
# Store my plots
HealthyPlots <- MyPlots


# Periodontitis Plots
PerioPlots <- (PerioPlots[[1]] | PerioPlots[[2]] | PerioPlots[[3]]) + plot_annotation(title = "Periodontitis")

# Healthy Plots
HealthyPlots <- (HealthyPlots[[1]] | HealthyPlots[[2]] | HealthyPlots[[3]]) + plot_annotation(title = "Healthy")

# Combine Healthy and Periodontitis
DA_Smokers <- (HealthyPlots / PerioPlots) +  plot_layout(guides = 'collect') + plot_annotation(title = "Smokers") & theme(legend.position = 'right')
DA_Smokers

ggsave(plot = DA_Smokers, 
       device = cairo_pdf,
       width = 11, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/DA_Smokers.pdf")



DA_Overall | DA_Smokers

DA_Overall








# ================================================================================================================================


set.seed(131219)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")
pseq_genus <- subset_samples(pseq_genus, Smoking == "No") # Subletting only smoker individuals

MySite <- unique(pseq_genus@sam_data$Site)

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


res$Group = "Periodontitis - Non Smokers"
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
                         #norm = "TSS",
                         method = "welch.test",
                         #p_adjust = "fdr",
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


res$Group = "Healthy - Non Smokers"
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
Periodontitis <- subset(results, Group == "Periodontitis - Non Smokers" & SigGenera == 1)
Healthy <- subset(results, Group == "Healthy - Non Smokers"  & SigGenera == 1)


# ========================================================#
#                                                         # 
#       ---- >    START PLOTTING HERE     < ----          #
#                                                         #   
# ========================================================#

# Define my theme
MyTheme <- theme_bw(12)+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Avenir", color = "black"),
        aspect.ratio = 1,
        legend.position = "bottom")

# Defining dimention of circle
MyLimits <- c(-1.5, 5.5)

# Defining my color palette
# this is a color-blindness friendly palette Ref. http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
#MyPalette <-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


MyPalette <- viridis(option = "E", 5)

# Setup datasets for cycle
Site <- c("Saliva", "Plaque", "Subgingival")
MyPlots <- list()

# ===============================
# Generate plot for Periodontitis
# ===============================
for (s in Site) {
  
  
  foo <- subset(Periodontitis, Site == s)
  
  # Compute number of agreement
  m <- ddply(foo, c("Genus", "Site"), summarize,
             Agreement = sum(!is.na(SigGenera)))
  foo <- merge(foo, m)
  
  
  
  # Get the name and the y position of each label
  label_data <- foo
  label_data <- subset(label_data[,c(1,7)], !duplicated(label_data$Genus))
  
  # Sorting data
  label_data <-label_data[order(-label_data$Agreement),]
  foo <-foo[order(-foo$Agreement),]
  
  # Creating custom id
  label_data$id <- seq(1, nrow(label_data))
  
  # get numer of bars
  number_of_bar <- nrow(label_data)
  
  # 
  angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #foo <- merge(foo, label_data)
  
  p <- ggplot(foo, aes(x=reorder(Genus, -Agreement), y=SigGenera))+
    geom_bar(stat = "identity", aes(fill = Model), width = .9)+
    scale_y_continuous(limits = MyLimits)+                # Circle dimension
    scale_fill_manual(values = MyPalette, drop = FALSE)+   # My color palette
    ggtitle("Periodontitis - Non Smokers", subtitle=s)+                                           # My title
    coord_polar()+
    # Add labels on top of each bar
    geom_text(data = label_data, aes(x=reorder(Genus, -Agreement), y=Agreement+0.02, label=Genus, hjust=hjust), 
              family = "Avenir", size=2.47, # 7Pt * Magic number 0.352777778
              angle = label_data$angle, 
              inherit.aes = FALSE) +
    MyTheme +
    theme_void()
  
  MyPlots[[paste0("Periodontitis-Non Smokers", s)]] <- p
  
  rm(foo, label_data, m, p)
  
}
rm(s)
# Store my plots
PerioPlots <- MyPlots


# reset MyPlots
MyPlots <- list()

# ===============================
# Generate plot for Healthy
# ===============================

for (s in Site) {
  
  foo <- subset(Healthy, Site == s)
  
  # Compute number of agreement
  m <- ddply(foo, c("Genus", "Site"), summarize,
             Agreement = sum(!is.na(SigGenera)))
  foo <- merge(foo, m)
  
  
  
  # Get the name and the y position of each label
  label_data <- foo
  label_data <- subset(label_data[,c(1,7)], !duplicated(label_data$Genus))
  
  # Sorting data
  label_data <-label_data[order(-label_data$Agreement),]
  foo <-foo[order(-foo$Agreement),]
  
  # Creating custom id
  label_data$id <- seq(1, nrow(label_data))
  
  # get numer of bars
  number_of_bar <- nrow(label_data)
  
  # 
  angle <- 90 - 360 * (label_data$id-0.5)/number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  #foo <- merge(foo, label_data)
  
  p <- ggplot(foo, aes(x=reorder(Genus, -Agreement), y=SigGenera))+
    geom_bar(stat = "identity", aes(fill = Model), width = .9)+
    scale_y_continuous(limits = MyLimits)+                # Circle dimension
    scale_fill_manual(values = MyPalette, drop = FALSE)+   # My color palette
    ggtitle("Healthy - Non Smokers", subtitle=s)+                                           # My title
    coord_polar()+
    # Add labels on top of each bar
    geom_text(data = label_data, aes(x=reorder(Genus, -Agreement), y=Agreement+0.02, label=Genus, hjust=hjust), 
              family = "Avenir", size=2.47, # 7Pt * Magic number 0.352777778
              angle = label_data$angle, 
              inherit.aes = FALSE) +
    MyTheme +
    theme_void()
  
  MyPlots[[paste0("Healthy-Non Smokers", s)]] <- p
  
  rm(foo, label_data, m, p)
  
}
rm(s)
# Store my plots
HealthyPlots <- MyPlots


# Periodontitis Plots
PerioPlots <- (PerioPlots[[1]] | PerioPlots[[2]] | PerioPlots[[3]])

# Healthy Plots
HealthyPlots <- (HealthyPlots[[1]] | HealthyPlots[[2]] | HealthyPlots[[3]])

# Combine Healthy and Periodontitis
DA_NonSmokers <- (HealthyPlots / PerioPlots) +  plot_layout(guides = 'collect') + plot_annotation(title = "Non Smokers") & theme(legend.position = 'right')
DA_NonSmokers

ggsave(plot = DA_NonSmokers, 
       device = cairo_pdf,
       width = 18, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/DA_NonSmokers.pdf")




























