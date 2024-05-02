library(phyloseq)
library(ggpubr)
library(ggplot2)
library(microbiome)
library(plyr)
library(reshape2)
library(microbiomeMarker)
library(extrafont)
library(showtext)
loadfonts(device = "pdf")

# Load my 16S data
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Setting my variables
Sites <- c("Saliva", "Plaque", "Subgingival")
res <- data.frame()

# # ================================
# #
# # Clycle for periodontitis
# 
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
# #
# # Clycle for Healthy
# 
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


# Asterisc if p ajusted <0.05
res$l <- ifelse(res$padj <0.05, "*", NA)


# Keep only what i need
ToKeep <- ddply(res, c("feature", "Site"), summarise,
                Is_significant = sum(!is.na(l)))
ToKeep <- subset(ToKeep, Is_significant >0)


# Mytheme
MyTheme <- theme_minimal(base_size = 11) + 
                  theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black", size = 7),
                  axis.text.x = element_text(angle = 35, hjust = 1, family = "Avenir"),
                  axis.text.y = element_text(family = "Avenir", size = 7),
                  text = element_text(family = "Avenir"),
                  aspect.ratio = 2/.25)



# reorder Site
res$Site <- factor(res$Site, levels = c("Saliva", "Plaque", "Subgingival"))
PerioHeatmap <-
 ggplot(data=subset(res, Periodontitis == "Periodontitis" & feature %in% ToKeep$feature), 
                          aes(y=feature, x=BioProject, fill=ef_diff_mean)) + 
                                geom_tile() + 
                                scale_x_discrete(drop=FALSE)+
                                scale_fill_gradient2(midpoint=0, high="#f2acda", mid="white",low="#6c9bd9", space ="Lab")+
                                geom_text(aes(label=l, family = "Avenir"), hjust = .5, ) +
                                facet_wrap(~Site, scales="free_x", drop = FALSE)+
                                labs(x="", y="", title = "Periodontitis")+
                                MyTheme


PerioHeatmap
ggsave(plot = PerioHeatmap,
       device = cairo_pdf,
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Heatmap_Periodontitis.pdf")









# reorder Site
res$Site <- factor(res$Site, levels = c("Saliva", "Plaque", "Subgingival"))
HealthyHeatmap <- ggplot(data=subset(res, Periodontitis == "Healthy"), 
                             aes(y=feature, x=BioProject, fill=ef_diff_mean)) + 
                      geom_tile() + 
                      scale_x_discrete(drop=FALSE)+
                      scale_fill_gradient2(midpoint=0, high="#f2acda", mid="white",low="#6c9bd9", space ="Lab")+
                      geom_text(aes(label=l, family = "Avenir"), hjust = .5, ) +
                      facet_wrap(~Site, scales="free_x", drop = FALSE)+
                      labs(x="", y="", title = "Healthy")+
                      MyTheme +
                      theme(aspect.ratio = 12/1,
                            axis.text.y = element_text(family = "Avenir", size = 5))


HealthyHeatmap
ggsave(plot = HealthyHeatmap,
       device = cairo_pdf,
       width = 9,
       filename = "AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Heatmap_Healthy.pdf")


  
