library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library(ggprism)

names(pdfFonts())

loadfonts(device = "pdf")



load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")


# define my smoking factors
psm_healthy@sam_data$Smoking <- ifelse(is.na(psm_healthy@sam_data$Smoking), 'Undefined', psm_healthy@sam_data$Smoking)
psm_healthy@sam_data$Smoking <- factor(psm_healthy@sam_data$Smoking, levels = c('No', 'Yes', 'Undefined'), labels = c('Non-smokers', 'Smokers', 'Undefined'))
psm_healthy@sam_data$Site <- factor(psm_healthy@sam_data$Site , levels = c('Saliva', 'Plaque', 'Subgingival'))
psm_healthy@sam_data$Sex <- factor(psm_healthy@sam_data$Sex , levels = c('female', 'male'), labels = c("F", "M"))


psm_perio@sam_data$Smoking <- ifelse(is.na(psm_perio@sam_data$Smoking), 'Undefined', psm_perio@sam_data$Smoking)
psm_perio@sam_data$Smoking <- factor(psm_perio@sam_data$Smoking, levels = c('No', 'Yes', 'Undefined'), labels = c('Non-smokers', 'Smokers', 'Undefined'))
psm_perio@sam_data$Site <- factor(psm_perio@sam_data$Site , levels = c('Saliva', 'Plaque', 'Subgingival'))
psm_perio@sam_data$Sex <- factor(psm_perio@sam_data$Sex , levels = c('female', 'male'), labels = c("F", "M"))



MyTheme <-  theme_classic(base_size = 10)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 2/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

set.seed(101)
p1 <- plot_richness(subset_samples(psm_perio, Smoking != "Undefined"), x="Sex", measures="Shannon", color = "Sex")+
  facet_grid(Site~Smoking)+
  geom_jitter(width = .15, alpha = I(.4))+
  geom_boxplot(width =.65, fill = NA, outlier.shape = NA)+
  labs(title = "Periodontitis")+
  coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme
p1$layers[1] <- NULL


p2 <- plot_richness(psm_perio, x="Sex", measures="Shannon", color = "Sex")+
  facet_grid(Site~.)+
  geom_jitter(width = .15, alpha = I(.4))+
  geom_boxplot(width =.65, fill = NA, outlier.shape = NA)+
  labs(title = "Periodontitis", subtitle = "Overall")+
  coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme 
p2$layers[1] <- NULL



p3 <- plot_richness(subset_samples(psm_healthy, Smoking != "Undefined"), x="Sex", measures="Shannon", color = "Sex")+
  facet_grid(Site~Smoking)+
  geom_jitter(width = .15, alpha = I(.4))+
  geom_boxplot(width =.65, fill = NA, outlier.shape = NA)+
  labs(title = "Healthy")+
  coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme 
p3$layers[1] <- NULL


p4 <- plot_richness(psm_healthy, x="Sex", measures="Shannon", color = "Sex")+
  facet_grid(Site~.)+
  geom_jitter(width = .15, alpha = I(.4))+
  geom_boxplot(width =.65, fill = NA, outlier.shape = NA)+
  labs(title = "Healthy")+
  coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme 
p4$layers[1] <- NULL


library(patchwork)

Alpha <- p1 + p3




ggsave(plot = Alpha, 
       device = cairo_pdf,
       width = 8.5, 
       height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Alpha_diversity.pdf")




# Chao 1 index

set.seed(101)
p1 <- plot_richness(subset_samples(psm_perio, Smoking != "Undefined"), x="Sex", measures="Chao1", color = "Sex")+
  facet_grid(Site~Smoking)+
  geom_jitter(width = .15, alpha = I(.4))+
  geom_boxplot(width =.65, fill = NA, outlier.shape = NA)+
  labs(title = "Periodontitis")+
  #coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme
p1$layers[1] <- NULL

p1

p3 <- plot_richness(subset_samples(psm_healthy, Smoking != "Undefined"), x="Sex", measures="Chao1", color = "Sex")+
  facet_grid(Site~Smoking)+
  geom_jitter(width = .15, alpha = I(.4))+
  geom_boxplot(width =.65, fill = NA, outlier.shape = NA)+
  labs(title = "Healthy")+
  #coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme 
p3$layers[1] <- NULL
p3








p1 <- plot_richness(subset_samples(psm_perio, Site == "Subgingival"), 
                    x="Sex", measures="Shannon", color = "Sex")+
  #facet_grid(Site~Smoking)+
  geom_boxplot(width =.5)+
  labs(title = "Periodontitis", subtitle = "Subgingival")+
  #coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme 


p2 <- plot_richness(subset_samples(psm_healthy, Site == "Subgingival"), 
                    x="Sex", measures="Shannon", color = "Sex")+
  #facet_grid(Site~Smoking)+
  geom_boxplot(width =.5)+
  labs(title = "Healthy", subtitle = "Subgingival")+
  #coord_fixed(ylim = c(0,4))+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme 

p1+p2
