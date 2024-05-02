library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library(patchwork)

loadfonts(device = "pdf")


# load my dataset
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# Defining factors
psm_perio@sam_data$Sex <- factor(psm_perio@sam_data$Sex, levels = c("female", "male"), labels = c("F", "M"))
psm_perio@sam_data$Site <- factor(psm_perio@sam_data$Site, levels = c("Saliva", "Plaque", "Subgingival"))

psm_healthy@sam_data$Sex <- factor(psm_healthy@sam_data$Sex, levels = c("female", "male"), labels = c("F", "M"))
psm_healthy@sam_data$Site <- factor(psm_healthy@sam_data$Site, levels = c("Saliva", "Plaque", "Subgingival"))


# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 3/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))


# Plot periodontitis
p1 <- plot_richness(psm_perio, x="Sex", measures="Observed", color = "Sex")+
          facet_grid(~Site)+
          geom_violin(aes(fill=Sex), alpha=.3)+
          geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=1.3, stroke =.4)+
          geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
                labs(title = "Periodontitis")+
                #coord_fixed(ylim = c(0,4))+
                scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
                scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
                MyTheme+
  stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family="Avenir")
  

# Remove dots
p1$layers[1] <- NULL
#p1


# Plot healthy
p0 <- plot_richness(psm_healthy, x="Sex", measures="Observed", color = "Sex")+
              facet_grid(~Site)+
              geom_violin(aes(fill=Sex), alpha=.3)+
              geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=1.3, stroke =.4)+
              geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
                    labs(title = "Healthy")+
                    #coord_fixed(ylim = c(0,4))+
                    scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
                    scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
                    MyTheme+
          stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family="Avenir")


# Remove dots
p0$layers[1] <- NULL
#p0

# Combining plot
Alpha <- p0+p1 + plot_layout(ncol = 1)

Alpha

ggsave(plot = Alpha, 
       device = cairo_pdf,
       width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Alpha_diversity.pdf")



