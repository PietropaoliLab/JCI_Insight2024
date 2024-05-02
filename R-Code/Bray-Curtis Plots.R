load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM PlyloSeq Object for analysys.RData")
library(phyloseq)
library(ggplot2)

set.seed(101)
ps_rarefied <- rarefy_even_depth(psm_ps)
ps_rarefied@sam_data$Sex <- factor(ps_rarefied@sam_data$Sex, levels = c("male", "female"), labels = c("Male", "Female"))

library(dplyr)
Smokers <- subset_samples(ps_rarefied, Smoking == "Yes")
Nonsmokers <- subset_samples(ps_rarefied, Smoking == "No")



MyTheme <- theme_bw(12)+
                theme(axis.text=element_text(size=10, color = "black"), panel.grid=element_blank(),
                strip.background = element_rect(colour=NA,fill="grey95"),
                aspect.ratio = 1,
                legend.title = element_blank())


p <- ordinate(ps_rarefied, "PCoA", "bray") %>% 
                plot_ordination(ps_rarefied, ., color = "Sex", title = "Bray-Curtis")+
                geom_point(aes(fill = Sex, shape = Site), size = 4)+
                facet_grid(Periodontitis ~ Smoking, scales = "free")+
                scale_fill_manual(values=c("#8ECEFD", "#F88B9D"))+
                scale_color_manual(values=c("#8ECEFD", "#F88B9D"))+
MyTheme

p
ggsave("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Bray-Curtis.eps", 
       plot = p, width = 10, height = 5.5, units = "in", dpi = 300)






# Make plots SMOKERS
MySite <- c("Plaque", "Saliva", "Subgingival")
MyPerio <- c("Healthy", "Periodontitis")
i = 1
Smokers_plot_list = list()
  for (perio in MyPerio) {  
    for (site in MySite) {
  foo <- subset_samples(Smokers,  Periodontitis == perio & Site == site)
 p =  ordinate(foo, "PCoA", "bray") %>% 
    plot_ordination(foo, ., color = "Sex", title = paste0(perio, "-", site)) +
    scale_fill_manual(values=c("#8ECEFD", "#F88B9D"))+
    scale_color_manual(values=c("#8ECEFD", "#F88B9D"))+
    MyTheme
  
  Smokers_plot_list[[i]] = p
  i = i+1
  rm(foo)
}
}

Smokers_plot_list


# Make plots NON-SMOKERS
MySite <- c("Plaque", "Saliva", "Subgingival")
MyPerio <- c("Healthy", "Periodontitis")
i = 1

NonSmokers_plot_list = list()
for (perio in MyPerio) {  
  for (site in MySite) {
    foo <- subset_samples(Nonsmokers,  Periodontitis == perio & Site == site)
    p =  ordinate(foo, "PCoA", "bray") %>% 
      plot_ordination(foo, ., color = "Sex", title = paste0(perio, "-", site)) +
      scale_fill_manual(values=c("#8ECEFD", "#F88B9D"))+
      scale_color_manual(values=c("#8ECEFD", "#F88B9D"))+
      MyTheme
    
    NonSmokers_plot_list[[i]] = p
    i = i+1
    rm(foo)
  }
}

NonSmokers_plot_list


library(ggpubr)
Psmokers <- ggarrange(Smokers_plot_list[[1]], Smokers_plot_list[[2]], Smokers_plot_list[[3]], Smokers_plot_list[[4]], 
                      Smokers_plot_list[[5]], Smokers_plot_list[[6]],
                      nrow=2, ncol = 3, common.legend = TRUE, legend="bottom")
Psmokers <- annotate_figure(Psmokers, fig.lab = "Smoking", fig.lab.size = 16)


Pnonsmokers <- ggarrange(NonSmokers_plot_list[[1]], NonSmokers_plot_list[[2]], NonSmokers_plot_list[[3]], NonSmokers_plot_list[[4]], 
                      NonSmokers_plot_list[[5]], NonSmokers_plot_list[[6]],
                      nrow=2, ncol = 3, common.legend = TRUE, legend="bottom")
Pnonsmokers <- annotate_figure(Pnonsmokers, fig.lab = "Non smoking", fig.lab.size = 16)

ggsave("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Bray-Curtis Smokers.eps", 
       plot = Psmokers, width = 10, height = 5.5, units = "in", dpi = 300)

ggsave("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Bray-Curtis NON Smokers.eps", 
       plot = Pnonsmokers, width = 10, height = 5.5, units = "in", dpi = 300)

rm(list = ls())


