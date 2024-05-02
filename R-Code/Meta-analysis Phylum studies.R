library(tidyverse)
library(SIAMCAT)
library(ggplot2)
library(phyloseq)
library(microbiome)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

MyPhylum <- data.frame(tax_table(psm_perio)[,"Phylum"])
MyPhylum$NewRownames <- paste0("OTU_", sprintf('%0.3d', 1:nrow(MyPhylum)), "_Phylum_", MyPhylum$Phylum)

taxa_names(psm_perio) <- MyPhylum$NewRownames
taxa_names(psm_perio)
Sites <- c("Saliva", "Plaque", "Subgingival")
assoc.list <- list()
k=1

for (s in Sites) {
  foo <- subset_samples(psm_perio, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  # TUTORIAL
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_meta.html
  for (d in datasets){
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    ps <- tax_glom(ps, taxrank = "Phylum")
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
    temp$PhylumRowname <- rownames(temp)
    temp$Phylum <- sc.obj@phyloseq@tax_table[,2]
    temp$Site = s
    
    assoc.list[[k]] <- temp %>% 
      select(Phylum, fc, auc, p.adj, Site) %>% 
      mutate(BioProject=d) # Change here
    k = k+1
  }
  
}
rm(temp, datasets, ps, sc.obj, d, label, s, k, foo, Sites)


# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(Phylum!='unclassified')

table(df.assoc$BioProject, df.assoc$Site)
table(psm_perio@sam_data$BioProject, psm_perio@sam_data$Site)


Phylum.of.interest <- df.assoc %>% 
  group_by(Phylum) %>% 
  summarise(m=mean(auc), n.filt=any(auc < 0.30 | auc > 0.70), 
            .groups='keep') %>% 
  filter(n.filt) %>% 
  arrange(m)


df_perio <- df.assoc %>% 
  # take only genera of interest
  filter(Phylum %in% Phylum.of.interest$Phylum) %>% 
  
  # convert to factor to enforce an ordering by mean AUC
  mutate(Phylum=factor(Phylum, levels = rev(Phylum.of.interest$Phylum))) %>% 
  
  # convert to factor to enforce ordering again
  mutate(BioProject=factor(BioProject)) %>% 
  
  # annotate the cells in the heatmap with stars
  mutate(l=case_when(p.adj < 0.05~'*', TRUE~'')) 





# Healthy

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

MyPhylum <- data.frame(tax_table(psm_healthy)[,"Phylum"])
MyPhylum$NewRownames <- paste0("OTU_", sprintf('%0.3d', 1:nrow(MyPhylum)), "_Phylum_", MyPhylum$Phylum)

taxa_names(psm_healthy) <- MyPhylum$NewRownames
taxa_names(psm_healthy)
Sites <- c("Saliva", "Plaque", "Subgingival")

k=1
for (s in Sites) {
  foo <- subset_samples(psm_healthy, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  # TUTORIAL
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_meta.html
  for (d in datasets){
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    ps <- tax_glom(ps, taxrank = "Phylum")
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
    temp$PhylumRowname <- rownames(temp)
    temp$Phylum <- sc.obj@phyloseq@tax_table[,2]
    temp$Site = s
    
    assoc.list[[k]] <- temp %>% 
      select(Phylum, fc, auc, p.adj, Site) %>% 
      mutate(BioProject=d) # Change here
    k = k+1
  }
  
}



rm(temp, datasets, ps, sc.obj, d, label, s, k, foo, Sites)


# combine all associations
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(Phylum!='unclassified')

table(df.assoc$BioProject, df.assoc$Site)
table(psm_perio@sam_data$BioProject, psm_perio@sam_data$Site)


Phylum.of.interest <- df.assoc %>% 
                        group_by(Phylum) %>% 
                            summarise(m=mean(auc), n.filt=any(auc < 0.30 | auc > 0.70), 
                            .groups='keep') %>% 
                            filter(n.filt) %>% 
                            arrange(m)


df_healthy <- df.assoc %>% 
  # take only genera of interest
  filter(Phylum %in% Phylum.of.interest$Phylum) %>% 
  
  # convert to factor to enforce an ordering by mean AUC
  mutate(Phylum=factor(Phylum, levels = rev(Phylum.of.interest$Phylum))) %>% 
  
  # convert to factor to enforce ordering again
  mutate(BioProject=factor(BioProject)) %>% 
  
  # annotate the cells in the heatmap with stars
  mutate(l=case_when(p.adj < 0.05~'*', TRUE~'')) 



p1 <- ggplot(data=df_perio, aes(y=Phylum, x=BioProject, fill=fc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
                       limits=c(-2.5, 3.5), name='Generalized\nfold change') + 
  theme_minimal() + 
  geom_text(aes(label=l), hjust = .5) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 35, hjust = 1)) + 
  facet_wrap(~Site, scales="free_x")+
  labs(x="", y="", title = "Periodontitis")


p2 <- ggplot(data=df_healthy, aes(y=Phylum, x=BioProject, fill=fc)) + 
  geom_tile() + 
  scale_fill_gradient2(low = '#3B6FB6', high='#D41645', mid = 'white', 
                       limits=c(-2.5, 3.5), name='Generalized\nfold change') + 
  theme_minimal() + 
  geom_text(aes(label=l), hjust = .5) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 35, hjust = 1)) + 
  facet_wrap(~Site, scales="free_x")+
  labs(x="", y="", title = "Healthy")





library(ggpubr)
ggarrange(p1, p2, common.legend = TRUE)



