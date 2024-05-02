library(phyloseq) 
library(microbiome)
library(caret)
library(MLeval)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(microbiomeMarker)
library(extrafont)
library(showtext)
loadfonts(device = "pdf")


# OSX
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")


######################################################
#
#         WITHIN STUDIES SAMPLE SIZE ANALYSIS
#
######################################################

healthy <- psm_healthy@sam_data
perio <- psm_perio@sam_data

heathy_perio <- rbind(healthy, perio)

Freq <- ddply(heathy_perio, c('BioProject','Site', 'Periodontitis'), summarise,
              Freqs = sum(!is.na(BioSample)))

Freq$Label <- paste0(Freq$BioProject, " - ", Freq$Periodontitis)
Freq$LongLabel <- paste0(Freq$BioProject, " - ", Freq$Periodontitis, " - ", Freq$Site)

# sorting
attach(Freq)
Freq <- Freq[order(Site, Periodontitis, Freqs),]
detach(Freq)

Freq$Ordering <- seq(1:nrow(Freq))

rm(heathy_perio, perio, healthy)



######################################################
#
#    EVALUATING GENERA DIFFERENCES WITHIN STUDIES
#           USING welch.test AND FDR
#
######################################################

# Setting my variables
Sites <- c("Saliva", "Plaque", "Subgingival")
res <- data.frame()

# # ================================
# # Clycle for periodontitis
# # ================================
# Within Sites
for (s in Sites) {
  foo <- subset_samples(psm_perio, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  
  # Within dataset
  for (d in datasets){
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Print some info on screen
    cat(paste0("Analyzing: ", d, "-", s, "\n"))
    
    # Operate my analysis
    tmp <- run_simple_stat(ps,
                           taxa_rank = "Genus",
                           transform = "log10p",
                           norm = "TSS",
                           method = "welch.test",
                           # p_adjust = "fdr",
                           pvalue_cutoff = 1, # Incorporate all Phyla
                           group = "Sex")
    # My results dataset
    tmp_res <- tmp@marker_table
    
    # Adding FDR to p-value
    tmp_res$padj <- p.adjust(tmp_res$pvalue, method = "fdr")
    
    # Adding type of test used
    tmp_res$Stat_method <- tmp@diff_method
    
    # Adding groups
    tmp_res$BioProject <- d
    tmp_res$Site <- s
    tmp_res$Periodontitis <- "Periodontitis"
    
    # Storing results
    res <- rbind(res, tmp_res)
    rm(tmp_res, tmp)
    
  }
}
rm(foo, ps, d, s)



# # ================================
# # Clycle for Healthy
# # ================================
# Within Sites
for (s in Sites) {
  foo <- subset_samples(psm_healthy, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  
  # Within dataset
  for (d in datasets){
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Print some info on screen
    cat(paste0("Analyzing: ", d, "-", s, "\n"))
    
    # Operate my analysis
    tmp <- run_simple_stat(ps,
                           taxa_rank = "Genus",
                           transform = "log10p",
                           norm = "TSS",
                           method = "welch.test",
                           # p_adjust = "fdr",
                           pvalue_cutoff = 1, # Incorporate all Phyla
                           group = "Sex")
    # My results dataset
    tmp_res <- tmp@marker_table
    
    # Adding FDR to p-value
    tmp_res$padj <- p.adjust(tmp_res$pvalue, method = "fdr")
    
    # Adding type of test used
    tmp_res$Stat_method <- tmp@diff_method
    
    # Adding groups
    tmp_res$BioProject <- d
    tmp_res$Site <- s
    tmp_res$Periodontitis <- "Healthy"
    
    # Storing results
    res <- rbind(res, tmp_res)
    rm(tmp_res, tmp)
    
  }
}
rm(foo, ps, d, s, Sites, datasets) # some cleaning



# Operate on welch.test dataset
res$LongLabel <- paste0(res$BioProject, " - ", res$Periodontitis, " - ", res$Site)
res$Label <- paste0(res$BioProject, " - ", res$Periodontitis)
res <- merge(res, Freq[,c('LongLabel', 'Ordering')], all.y = TRUE)

# Compute Significant Genera
welch <- ddply(res, c("LongLabel", "Label", "Ordering", "enrich_group", "Site"), summarise,
               N = sum(padj <0.05))

welch$Sex <- ifelse(welch$enrich_group == "female", "F", "M")


welch$Site <- factor(welch$Site, levels = c("Saliva", "Plaque", "Subgingival"))


rm(res)

MyDataset <- merge(Freq, welch, by ="LongLabel")


library(ggpubr)

ggscatter(MyDataset, x = "N", y = "Freqs", color = "Sex",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  labs(x="Significant enriched genera (N)", y="Study sample size (N)")+
  stat_cor(method = "pearson")  # Add correlation coefficient












