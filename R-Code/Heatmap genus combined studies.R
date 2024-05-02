library(tidyverse)
library(SIAMCAT)
library(ggplot2)
library(phyloseq)
library(microbiome)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")


Sites <- c("Saliva", "Plaque", "Subgingival")
assoc.list <- list()
k=1

for (s in Sites) {
  
  ps <- subset_samples(psm_perio, Site == s)
  
  # Get bioproject with specific site 
  #datasets <- unique(foo@sam_data$BioProject)
  
  # TUTORIAL
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_meta.html
  #for (d in datasets){
  #  ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Convert in relative abundance
    ps <- transform(ps, transform = "compositional")
    
    
    label <- create.label(meta=sample_data(ps),
                          label = "Sex",
                          case = "female")
    
    # run the constructor function
    sc.obj <- siamcat(phyloseq=ps, label=label)
    
    
    # test for associations
    sc.obj <- check.associations(sc.obj,  test='wilcoxon', feature.type = 'original')
    
    # extract the associations and save them in the assoc.list
    temp <- associations(sc.obj)
    temp$genus <- rownames(temp)
    temp$Site = s
    
    assoc.list[[k]] <- temp %>% 
      select(genus, fc, auc, p.adj) %>% 
      mutate(Site=s) # Change here
    k = k+1
  }
  

rm(temp, datasets, ps, sc.obj, d, label, s, k, foo, Sites)


# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(genus!='unclassified')

table(df.assoc$BioProject, df.assoc$Site)
table(psm_perio@sam_data$BioProject, psm_perio@sam_data$Site)


genera.of.interest <- df.assoc %>% 
  group_by(genus) %>% 
  summarise(m=mean(auc), n.filt=any(auc < 0.30 | auc > 0.70), 
            .groups='keep') %>% 
  filter(n.filt) %>% 
  arrange(m)


df <- df.assoc %>% 
  # take only genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  
  # convert to factor to enforce an ordering by mean AUC
  mutate(genus=factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  
  # convert to factor to enforce ordering again
  mutate(Site=factor(Site)) %>% 
  
  # annotate the cells in the heatmap with stars
  mutate(l=case_when(p.adj < 0.05~'*', TRUE~'')) 

# Make labs nicer
df$genus <- gsub(".*Genus_", "", df$genus)

p1 <- ggplot(data=df, aes(y=genus, x=Site, fill=fc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
                       limits=c(-2.5, 3.5), name='Generalized\nfold change') + 
  theme_minimal() + 
  geom_text(aes(label=l), hjust = .5) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 35, hjust = 1)) + 
  labs(x="", y="", title = "Periodontitis")

p1

# ############################################################
# HEALTHY
# ############################################################
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")


Sites <- c("Saliva", "Plaque", "Subgingival")
assoc.list <- list()
k=1

for (s in Sites) {
  
  ps <- subset_samples(psm_healthy, Site == s)
  
  # Get bioproject with specific site 
  #datasets <- unique(foo@sam_data$BioProject)
  
  # TUTORIAL
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_meta.html
  #for (d in datasets){
  #  ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
  
  # Convert in relative abundance
  ps <- transform(ps, transform = "compositional")
  
  
  label <- create.label(meta=sample_data(ps),
                        label = "Sex",
                        case = "female")
  
  # run the constructor function
  sc.obj <- siamcat(phyloseq=ps, label=label)
  
  
  # test for associations
  sc.obj <- check.associations(sc.obj,  test='wilcoxon', feature.type = 'original')
  
  # extract the associations and save them in the assoc.list
  temp <- associations(sc.obj)
  temp$genus <- rownames(temp)
  temp$Site = s
  
  assoc.list[[k]] <- temp %>% 
    select(genus, fc, auc, p.adj) %>% 
    mutate(Site=s) # Change here
  k = k+1
}


rm(temp, datasets, ps, sc.obj, d, label, s, k, foo, Sites)


# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(genus!='unclassified')

table(df.assoc$BioProject, df.assoc$Site)
table(psm_perio@sam_data$BioProject, psm_perio@sam_data$Site)


genera.of.interest <- df.assoc %>% 
  group_by(genus) %>% 
  summarise(m=mean(auc), n.filt=any(auc < 0.30 | auc > 0.70), 
            .groups='keep') %>% 
  filter(n.filt) %>% 
  arrange(m)


df <- df.assoc %>% 
  # take only genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  
  # convert to factor to enforce an ordering by mean AUC
  mutate(genus=factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  
  # convert to factor to enforce ordering again
  mutate(Site=factor(Site)) %>% 
  
  # annotate the cells in the heatmap with stars
  mutate(l=case_when(p.adj < 0.05~'*', TRUE~'')) 

# Make labs nicer
df$genus <- gsub(".*Genus_", "", df$genus)


p2 <- ggplot(data=df, aes(y=genus, x=Site, fill=fc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
                       limits=c(-2.5, 3.5), name='Generalized\nfold change') + 
  theme_minimal() + 
  geom_text(aes(label=l), hjust = .5) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 35, hjust = 1)) + 
  labs(x="", y="", title = "Healthy")

p2

library(ggpubr)
ggarrange(p1, p2, common.legend = TRUE)



