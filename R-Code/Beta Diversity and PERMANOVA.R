library("microbial")
library(phyloseq)
library(microbiome)

# preRef(ref_db = "silva", path=".")

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")




psm_healthy@sam_data$Smoking <- ifelse(is.na(psm_healthy@sam_data$Smoking), 'Undefined', psm_healthy@sam_data$Smoking)
psm_healthy@sam_data$Smoking <- factor(psm_healthy@sam_data$Smoking, levels = c('No', 'Yes', 'Undefined'), labels = c('Non-smokers', 'Smokers', 'Undefined'))

psm_perio@sam_data$Smoking <- ifelse(is.na(psm_perio@sam_data$Smoking), 'Undefined', psm_perio@sam_data$Smoking)
psm_perio@sam_data$Smoking <- factor(psm_perio@sam_data$Smoking, levels = c('No', 'Yes', 'Undefined'), labels = c('Non-smokers', 'Smokers', 'Undefined'))


MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))

# Define my smoking habits for cycle
MySmoking <- unique(psm_healthy@sam_data$Smoking)

# Define my resuls object
MyPlots <- list()

# define my index for plots array
k = 1

# Beta Diversity
for (s in MySmoking) {
  
  # Subsetting data of interest
  tmp <- subset_samples(psm_healthy, Smoking == s)
  
foo <- plotbeta(tmp, group="Site", 
                distance = "bray", 
                ellipse = FALSE, 
                color = c("#7fc97f", "#1f78b4", "#d01c8b"))+
                MyTheme + 
                labs(title = "Healthy", color = "Site", subtitle = s) + theme(aspect.ratio = 1)
# add alpha               
foo$layers[[1]]$aes_params$alpha <- .3

MyPlots[[k]] <- foo
k=k+1
rm(tmp, foo)

}




# Beta Diversity
for (s in MySmoking) {
  
  # Subsetting data of interest
  tmp <- subset_samples(psm_perio, Smoking == s)
  
  foo <- plotbeta(tmp, group="Site", 
                  distance = "bray", 
                  ellipse = FALSE, 
                  color = c("#7fc97f", "#1f78b4", "#d01c8b"))+
    MyTheme + 
    labs(title = "Periodontitis", color = "Site", subtitle = s) + theme(aspect.ratio = 1)
  # add alpha               
  foo$layers[[1]]$aes_params$alpha <- .3
  
  MyPlots[[k]] <- foo
  k=k+1
  rm(tmp, foo)
  
}




# Plotting results
library(patchwork)
MyPlots[[1]] + MyPlots[[2]] + MyPlots[[3]] + 
MyPlots[[4]] + MyPlots[[5]] + MyPlots[[6]] + plot_layout(nrow = 2, guides = 'collect')





tmp <- subset_samples(psm_perio, Site == "Subgingival")
p1 <- plotbeta(tmp, group="Sex", 
         distance = "bray", 
         ellipse = FALSE, 
         color = c("#7fc97f", "#1f78b4", "#d01c8b"))+
  MyTheme + 
  labs(title = "Periodontitis", color = "Sex", subtitle = "Subgingival") + theme(aspect.ratio = 1)


tmp <- subset_samples(psm_healthy, Site == "Subgingival")
p2 <- plotbeta(tmp, group="Sex", 
               distance = "bray", 
               ellipse = FALSE, 
               color = c("#7fc97f", "#1f78b4", "#d01c8b"))+
  MyTheme + 
  labs(title = "Healthy", color = "Sex", subtitle = "Subgingival") + theme(aspect.ratio = 1)

library(patchwork)
p1+p2







# PERMANOVA
# https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
library(vegan)
set.seed(101)


# Periodontitis

# Calculate bray curtis distance matrix
tmp <- subset_samples(psm_perio, Smoking == 'Non-smoker')
erie_bray <- phyloseq::distance(tmp, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psm_perio))

# Adonis test
adonis2(erie_bray ~ Site, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(erie_bray, sampledf$Site)
permutest(beta)




