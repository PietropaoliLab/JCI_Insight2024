library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library(patchwork)

loadfonts(device = "pdf")


# load my dataset
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# Defining factors
psm_perio@sam_data$Sex <- factor(psm_perio@sam_data$Sex, levels = c("female", "male"), labels = c("F", "M"))
psm_perio@sam_data$Site <- factor(psm_perio@sam_data$Site, levels = c("Saliva", "Plaque", "Subgingival"))
psm_perio <- subset_samples(psm_perio, !is.na(Smoking))
psm_perio@sam_data$Smoking <- factor(psm_perio@sam_data$Smoking, levels = c("Yes", "No"), labels = c("Smokers", "Non smokers"))

psm_healthy@sam_data$Sex <- factor(psm_healthy@sam_data$Sex, levels = c("female", "male"), labels = c("F", "M"))
psm_healthy@sam_data$Site <- factor(psm_healthy@sam_data$Site, levels = c("Saliva", "Plaque", "Subgingival"))
psm_healthy <- subset_samples(psm_healthy, !is.na(Smoking))
psm_healthy@sam_data$Smoking <- factor(psm_healthy@sam_data$Smoking, levels = c("Yes", "No"), labels = c("Smokers", "Non smokers"))


# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 2.5/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))


# Plot periodontitis
p1 <- plot_richness(psm_perio, x="Sex", measures="Observed", color = "Sex")+
                    facet_grid(Smoking~Site)+
                    geom_violin(aes(fill=Sex), alpha=.3)+
                    geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=2.4)+
                    geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
                    labs(title = "Periodontitis", y = "Observed AVS diversity")+
                    scale_y_continuous(breaks = c(20, 40, 60, 80))+
                    coord_fixed(ylim = c(10,90))+
                    scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
                    scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
                    MyTheme+
                    stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family="Avenir")
                  

# Remove dots
p1$layers[1] <- NULL
p1


# Plot healthy
p0 <- plot_richness(psm_healthy, x="Sex", measures="Observed", color = "Sex")+
                    facet_grid(Smoking~Site)+
                    geom_violin(aes(fill=Sex), alpha=.3)+
                    geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=2.4)+
                    geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
                    labs(title = "Healthy",  y = "Observed AVS diversity")+
                    scale_y_continuous(breaks = c(20, 40, 60, 80))+
                    coord_fixed(ylim = c(10,90))+
                    scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
                    scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
                    MyTheme+
                    stat_compare_means(comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family="Avenir")


# Remove dots
p0$layers[1] <- NULL
#p0

# Combining plot
Alpha <- (p0 | p1)

Alpha


ggsave(plot = Alpha, 
       device = cairo_pdf,
       width = 8.5, 
       height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Alpha_diversity_Smoking.pdf")









H <- p0$data
P <- p1$data

library(plyr)
library(patchwork)
library(ggpubr)

H$Log10Value <- log10(H$value + 1)
P$Log10Value <- log10(P$value + 1)


Hmu <- ddply(H, c("Sex", "Site", "Smoking"), summarise,
             grp.mean=mean(Log10Value))

Pmu <- ddply(P, c("Sex", "Site", "Smoking"), summarise,
             grp.mean=mean(Log10Value))



Ph <- ggplot(H, aes(x=Log10Value, color=Smoking, fill = Smoking)) +
              geom_density(alpha=0.25)+
              geom_vline(data=Hmu, aes(xintercept=grp.mean, color=Smoking),linetype="dashed", show.legend = FALSE)+
              scale_color_manual(values=c("#E69F00", "#56B4E9"))+
              scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
              labs(x="Observed ASV diversity", y="Density", title = "Healthy")+
              ylim(c(0, 12))+
              xlim(c(1, 2))+
              facet_grid(Sex~Site)
  

Pp <- ggplot(P, aes(x=Log10Value, color=Smoking, fill = Smoking)) +
                geom_density(alpha=0.25)+
                geom_vline(data=Pmu, aes(xintercept=grp.mean, color=Smoking),linetype="dashed", show.legend = FALSE)+
                scale_color_manual(values=c("#E69F00", "#56B4E9"))+
                scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
                ylim(c(0, 12))+
                xlim(c(1, 2))+
                labs(x="Observed ASV diversity", y="Density", title = "Periodontitis")+
                facet_grid(Sex~Site)


MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))



Density <- (Ph|Pp) + plot_layout(guides = "collect") & 
              MyTheme

Density
ggsave(plot = Density, 
       device = cairo_pdf,
       width = 7, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Density smokers.pdf")


compare_means(Log10Value ~ Smoking, group.by = c("Sex", "Site"), data = H)
compare_means(Log10Value ~ Smoking, group.by = c("Sex", "Site"), data = P)



