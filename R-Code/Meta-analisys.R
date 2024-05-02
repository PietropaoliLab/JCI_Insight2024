library(tidyverse)
library(SIAMCAT)
library(ggplot2)
library(phyloseq)
library(microbiome)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

datasets <- unique(psm_perio@sam_data$BioProject)

# TUTORIAL
# https://www.bioconductor.org/packages/devel/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_meta.html
assoc.list <- list()
for (d in datasets){
  ps <- subset_samples(psm_perio, BioProject == d)
  
  # Convert in relative abundance
  ps <- transform(ps, transform = "compositional")
  
  
  label <- create.label(meta=sample_data(ps),
                        label = "Sex",
                        case = "female")
  
  # run the constructor function
  sc.obj <- siamcat(phyloseq=ps, label=label)
  
  
  # test for associations
  sc.obj <- check.associations(sc.obj, feature.type = 'original')
  
  # extract the associations and save them in the assoc.list
  temp <- associations(sc.obj)
  temp$genus <- rownames(temp)
  assoc.list[[d]] <- temp %>% 
    select(genus, fc, auc, p.adj) %>% 
    mutate(BioProject=d)
}

# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(genus!='unclassified')




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
  mutate(BioProject=factor(BioProject, levels = datasets)) %>% 
  # annotate the cells in the heatmap with stars
  mutate(l=case_when(p.adj < 0.05~'*', TRUE~'')) 


  ggplot(data=df, aes(y=genus, x=BioProject, fill=fc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
                       limits=c(-2.7, 2.7), name='Generalized\nfold change') + 
  theme_minimal() + 
  geom_text(aes(label=l)) +
  theme(panel.grid = element_blank()) + 
  xlab('') + ylab('') +
  theme(axis.text = element_text(size=6))




