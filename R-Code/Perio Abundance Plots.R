# Load PSM healthy dataset
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")


library(microbiome)
library(phyloseq)
library(hrbrthemes)
library(ggpubr)

psm_perio@sam_data$Sex <- factor(psm_perio@sam_data$Sex, levels = c("female","male"), labels = c("F", "M"))

Plaque <- subset_samples(psm_perio, Site == "Plaque")
Saliva <- subset_samples(psm_perio, Site == "Saliva")
Subgingival <- subset_samples(psm_perio, Site == "Subgingival")

MyTheme <- theme_linedraw(12)+
                theme(axis.text = element_text(color = "black"),
                      axis.ticks.x=element_blank(),
                      panel.grid = element_blank(),
                      panel.border = element_blank())

# Averaged by group
# Use all samples 
# get relative abudance
Plaque.rel <- microbiome::transform(Plaque, "compositional")
Plaque.fam.rel <-aggregate_rare(Plaque.rel, level = "Phylum", detection = 0.005, prevalence = 0.5)

p <- plot_composition(Plaque.fam.rel,
                       average_by = "Sex") + 
  guides(fill = guide_legend(ncol = 1))+
  scale_fill_brewer("Phylum", palette = "Paired")+
  labs(x = "", 
       y = "Relative abundance",
       title = "Plaque")+
  scale_y_continuous(label = scales::percent) + 
  MyTheme
print(p)



  
  
# Averaged by group
# Use all samples 
# get relative abudance
Saliva.rel <- microbiome::transform(Saliva, "compositional")
Saliva.fam.rel <-aggregate_rare(Saliva.rel, level = "Phylum", detection = 0.005, prevalence = 0.5)

p1 <- plot_composition(Saliva.fam.rel,
                      average_by = "Sex") + 
  guides(fill = guide_legend(ncol = 1))+
  scale_fill_brewer("Phylum", palette = "Paired")+
  labs(x = "", 
       y = "Relative abundance",
       title = "Saliva")+
  scale_y_continuous(label = scales::percent) + 
  MyTheme
print(p1)





# Averaged by group
# Use all samples 
# get relative abudance
Subgingival.rel <- microbiome::transform(Subgingival, "compositional")
Subgingival.fam.rel <-aggregate_rare(Subgingival.rel, level = "Phylum", detection = 0.005, prevalence = 0.5)

p2 <- plot_composition(Subgingival.fam.rel,
                       average_by = "Sex") + 
  guides(fill = guide_legend(ncol = 1))+
  scale_fill_brewer("Phylum", palette = "Paired")+
  labs(x = "", 
       y = "Relative abundance",
       title = "Subgingival")+
  scale_y_continuous(label = scales::percent) + 
MyTheme
print(p2)


ggarrange(p, p1, p2, nrow = 1)
