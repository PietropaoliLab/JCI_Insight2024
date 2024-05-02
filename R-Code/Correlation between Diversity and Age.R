library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library(ggprism)

loadfonts(device = "pdf")

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

tmp <- merge_phyloseq(psm_perio, psm_healthy)

tmp@sam_data$Age <- as.numeric(tmp@sam_data$Age)
quantile(tmp@sam_data$Age, probs = c(0, .33, .66, 1))

tmp@sam_data$AgeCat <- cut(tmp@sam_data$Age, breaks = c(0, 40, 55, Inf), labels = c("<40", "40-54", "≥55"))


# define my smoking factors
tmp@sam_data$Smoking <- factor(tmp@sam_data$Smoking, levels = c('No', 'Yes'), labels = c('Non-smokers', 'Smokers'))



MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

set.seed(101)
p1 <- plot_richness(tmp, x="AgeCat", measures="Observed", color = "Sex")

# My Shannon Index
Diversity <- p1$data



Diversity$Site <- factor(Diversity$Site, levels = c("Saliva", "Plaque", "Subgingival"))

Diversity_plot <- ggplot(Diversity, aes(x = AgeYrs, y = value))+
                            geom_point(color = "gray90", shape = 19, alpha = .8, size = 2.2)+
                            geom_smooth(aes(color = Sex), method = "lm", se = FALSE, show.legend = FALSE, fullrange=TRUE)+
                            scale_color_manual(values = c("#f2acda", "#6c9bd9"))+
                            facet_grid(Periodontitis ~ Site)+
                            stat_cor(aes(color = Sex), label.y = c(7, 17), size = 4, family = "Avenir") +
                            labs(x="Age", y="Observed ASV Diversity")+
                            MyTheme+
                            coord_cartesian(expand = FALSE)



ggsave(plot = Diversity_plot, 
       device = cairo_pdf,
       width = 8, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Correlation Age vs Diversity.pdf")


