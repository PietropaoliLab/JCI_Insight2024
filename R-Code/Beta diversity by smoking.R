library("microbial")
library(phyloseq)
library(microbiome)
library(phyloseq)
library(reshape2)
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
psm_perio@sam_data$Smoking <- factor(psm_perio@sam_data$Smoking, levels = c("Yes", "No"), labels = c("Smokers", "Non smokers"))


psm_healthy@sam_data$Sex <- factor(psm_healthy@sam_data$Sex, levels = c("female", "male"), labels = c("F", "M"))
psm_healthy@sam_data$Site <- factor(psm_healthy@sam_data$Site, levels = c("Saliva", "Plaque", "Subgingival"))
psm_healthy@sam_data$Smoking <- factor(psm_healthy@sam_data$Smoking, levels = c("Yes", "No"), labels = c("Smokers", "Non smokers"))




# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 8),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

# Healthy Non smokers
HNS <- plotbeta(subset_samples(psm_healthy, Smoking == "Non smokers"), 
              group="Site", 
              distance = "bray", 
              ellipse = FALSE)+
              ylim (-0.5, 0.5)+
              xlim (-0.5, 0.5)+
              MyTheme
HNS

# Storing my x and y labels
Lab_x <- HNS$labels$x
Lab_y <- HNS$labels$y

# Storing my data
HNS <- HNS$data

# Put Axis labels into my dataset
HNS$Lab_x <- Lab_x
HNS$Lab_y <- Lab_y 
rm(Lab_x, Lab_y)

# Plotting beta diversity of healthy
KW <- melt(HNS[, c("Axis.1", "Axis.2", "Site")], id.vars = c("Site"))
pval <- paste0("P=", compare_means(value ~ Site, method = "kruskal.test", data = KW)$p.adj)
pval
PlotHNS <- ggplot(HNS, aes(x=Axis.1, y=Axis.2))+
                        geom_point(aes(color = Site), size = 3, alpha = .7)+
                        scale_color_manual(values = c("#ccece6", "#1c9099", "#006d2c"))+
                        labs(x=unique(HNS$Lab_x), y=unique(HNS$Lab_y), title = "Healthy - Non smokers")+
                        annotate(geom = "text", label = pval, x = .50, y = -.5, family = "Avenir", size = 2.8, hjust = 1)+
                        ylim (-0.5, 0.5)+
                        xlim (-0.5, 0.5)+
                        MyTheme

PlotHNS



# Healthy Smokers
HS <- plotbeta(subset_samples(psm_healthy, Smoking == "Smokers"), 
                group="Site", 
                distance = "bray", 
                ellipse = FALSE)+
  ylim (-0.5, 0.5)+
  xlim (-0.5, 0.5)+
  MyTheme
HS

# Storing my x and y labels
Lab_x <- HS$labels$x
Lab_y <- HS$labels$y

# Storing my data
HS <- HS$data

# Put Axis labels into my dataset
HS$Lab_x <- Lab_x
HS$Lab_y <- Lab_y 
rm(Lab_x, Lab_y)

# Plotting beta diversity of healthy
KW <- melt(HS[, c("Axis.1", "Axis.2", "Site")], id.vars = c("Site"))
pval <- paste0("P=", compare_means(value ~ Site, method = "kruskal.test", data = KW)$p.adj)
pval
PlotHS <- ggplot(HS, aes(x=Axis.1, y=Axis.2))+
              geom_point(aes(color = Site), size = 3, alpha = .7)+
              scale_color_manual(values = c("#ccece6", "#1c9099", "#006d2c"))+
              labs(x=unique(HS$Lab_x), y=unique(HS$Lab_y), title = "Healthy - Smokers")+
              annotate(geom = "text", label = pval, x = .50, y = -.5, family = "Avenir", size = 2.8, hjust = 1)+
              ylim (-0.5, 0.5)+
              xlim (-0.5, 0.5)+
              MyTheme

PlotHS



# Periodontitis - Non smokers
PNS <- plotbeta(subset_samples(psm_perio, Smoking == "Non smokers"), 
                group="Site", 
                distance = "bray", 
                ellipse = FALSE)+
                ylim (-0.5, 0.5)+
                xlim (-0.5, 0.5)+
                MyTheme
PNS

# Storing my x and y labels
Lab_x <- PNS$labels$x
Lab_y <- PNS$labels$y

# Storing my data
PNS <- PNS$data

# Put Axis labels into my dataset
PNS$Lab_x <- Lab_x
PNS$Lab_y <- Lab_y 
rm(Lab_x, Lab_y)

# Plotting beta diversity of healthy
KW <- melt(PNS[, c("Axis.1", "Axis.2", "Site")], id.vars = c("Site"))
pval <- paste0("P=", compare_means(value ~ Site, method = "kruskal.test", data = KW)$p.adj)
pval
PlotPNS <- ggplot(PNS, aes(x=Axis.1, y=Axis.2))+
                    geom_point(aes(color = Site), size = 3, alpha = .7)+
                    scale_color_manual(values = c("#fc9272", "#531C46", "#B9263A"))+
                    labs(x=unique(PNS$Lab_x), y=unique(PNS$Lab_y), title = "Periodontitis - Non smokers")+
                    annotate(geom = "text", label = pval, x = .50, y = -.5, family = "Avenir", size = 2.8, hjust = 1)+
                    ylim (-0.5, 0.5)+
                    xlim (-0.5, 0.5)+
                    MyTheme

PlotPNS



# Periodontitis - Smokers
PS <- plotbeta(subset_samples(psm_perio, Smoking == "Smokers"), 
                group="Site", 
                distance = "bray", 
                ellipse = FALSE)+
                ylim (-0.5, 0.5)+
                xlim (-0.5, 0.5)+
                MyTheme
PS

# Storing my x and y labels
Lab_x <- PS$labels$x
Lab_y <- PS$labels$y

# Storing my data
PS <- PS$data

# Put Axis labels into my dataset
PS$Lab_x <- Lab_x
PS$Lab_y <- Lab_y 
rm(Lab_x, Lab_y)

# Plotting beta diversity of healthy
KW <- melt(PS[, c("Axis.1", "Axis.2", "Site")], id.vars = c("Site"))
pval <- paste0("P=", compare_means(value ~ Site, method = "kruskal.test", data = KW)$p.adj)
pval
PlotPS <- ggplot(PS, aes(x=Axis.1, y=Axis.2))+
                geom_point(aes(color = Site), size = 3, alpha = .7)+
                scale_color_manual(values = c("#fc9272", "#531C46", "#B9263A"))+
                labs(x=unique(PS$Lab_x), y=unique(PS$Lab_y), title = "Periodontitis - Smokers")+
                annotate(geom = "text", label = pval, x = .50, y = -.5, family = "Avenir", size = 2.8, hjust = 1)+
                ylim (-0.5, 0.5)+
                xlim (-0.5, 0.5)+
                MyTheme

PlotPS




# ===================================== #

#     Compute P-values for Brays
#         Kruscal wallis
# ===================================== #



HealthyPalette <- c("#ccece6", "#1c9099", "#006d2c")
PerioPalette <- c("#fc9272", "#531C46", "#B9263A")

# Dim-1 Healthy Non smokers
pval <- paste0("P=", compare_means(Axis.1 ~ Site, method = "kruskal.test", data = HNS)$p.adj)
Dim1_HNS <- ggplot(HNS, aes(x = Site, y = Axis.1))+
                geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                scale_x_discrete(limits = rev(levels(HNS$Site)))+
                scale_fill_manual(values = HealthyPalette)+
                ylim(-0.5, 0.5)+
                coord_flip()+
                annotate(geom = "text", label = pval, x = 1, y= .45, family = "Avenir", size = 2.8)+
                labs(x="", y="Dim-1", title = "Healthy - Non smoker")+
                MyTheme+
                theme(aspect.ratio = 1/5,
                      legend.position = "none")
Dim1_HNS


# Dim-2 Healthy Non smokers
pval <- paste0("P=", compare_means(Axis.2 ~ Site, method = "kruskal.test", data = HNS)$p.adj)
Dim2_HNS <- ggplot(HNS, aes(x = Site, y = Axis.2))+
                  geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                  scale_x_discrete(limits = (levels(HNS$Site)))+
                  scale_fill_manual(values = HealthyPalette)+
                  ylim(-0.5, 0.5)+
                  annotate(geom = "text", label = pval, x = 2, y = .5, family = "Avenir", size = 2.8, hjust = 0.5)+
                  labs(x="", y="Dim-2", title = "Healthy - Non smoker")+
                  MyTheme+
                  theme(aspect.ratio = 5/1,
                        legend.position = "none")
Dim2_HNS





# Dim-1 Healthy Smokers
pval <- paste0("P=", compare_means(Axis.1 ~ Site, method = "kruskal.test", data = HS)$p.adj)
Dim1_HS <- ggplot(HS, aes(x = Site, y = Axis.1))+
                  geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                  scale_x_discrete(limits = rev(levels(HS$Site)))+
                  scale_fill_manual(values = HealthyPalette)+
                  ylim(-0.5, 0.5)+
                  coord_flip()+
                  annotate(geom = "text", label = pval, x = 1, y= .45, family = "Avenir", size = 2.8)+
                  labs(x="", y="Dim-1", title = "Healthy - Smoker")+
                  MyTheme+
                  theme(aspect.ratio = 1/5,
                        legend.position = "none")
Dim1_HS


# Dim-2 Healthy Smokers
pval <- paste0("P=", compare_means(Axis.2 ~ Site, method = "kruskal.test", data = HS)$p.adj)
Dim2_HS <- ggplot(HS, aes(x = Site, y = Axis.2))+
                      geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                      scale_x_discrete(limits = (levels(HS$Site)))+
                      scale_fill_manual(values = HealthyPalette)+
                      ylim(-0.5, 0.5)+
                      annotate(geom = "text", label = pval, x = 2, y = .5, family = "Avenir", size = 2.8, hjust = 0.5)+
                      labs(x="", y="Dim-2", title = "Healthy - Smoker")+
                      MyTheme+
                      theme(aspect.ratio = 5/1,
                            legend.position = "none")
Dim2_HS








# Dim-1 Perio Non smokers
pval <- paste0("P=", compare_means(Axis.1 ~ Site, method = "kruskal.test", data = PNS)$p.adj)
Dim1_PNS <- ggplot(PNS, aes(x = Site, y = Axis.1))+
                    geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                    scale_x_discrete(limits = rev(levels(PNS$Site)))+
                    scale_fill_manual(values = PerioPalette)+
                    ylim(-0.5, 0.5)+
                    coord_flip()+
                    annotate(geom = "text", label = pval, x = 1, y= .45, family = "Avenir", size = 2.8)+
                    labs(x="", y="Dim-1", title = "Periodontitis - Non smoker")+
                    MyTheme+
                    theme(aspect.ratio = 1/5,
                          legend.position = "none")
Dim1_PNS


# Dim-2 Periodontitis Non smokers
pval <- paste0("P=", compare_means(Axis.2 ~ Site, method = "kruskal.test", data = PNS)$p.adj)
Dim2_PNS <- ggplot(PNS, aes(x = Site, y = Axis.2))+
                    geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                    scale_x_discrete(limits = (levels(PNS$Site)))+
                    scale_fill_manual(values = PerioPalette)+
                    ylim(-0.5, 0.5)+
                    annotate(geom = "text", label = pval, x = 2, y = .5, family = "Avenir", size = 2.8, hjust = 0.5)+
                    labs(x="", y="Dim-2", title = "Periodontitis - Non smoker")+
                    MyTheme+
                    theme(aspect.ratio = 5/1,
                          legend.position = "none")
Dim2_PNS





# Dim-1 Periodontitis Smokers
pval <- paste0("P=", compare_means(Axis.1 ~ Site, method = "kruskal.test", data = PS)$p.adj)
Dim1_PS <- ggplot(PS, aes(x = Site, y = Axis.1))+
                  geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                  scale_x_discrete(limits = rev(levels(HS$Site)))+
                  scale_fill_manual(values = PerioPalette)+
                  ylim(-0.5, 0.5)+
                  coord_flip()+
                  annotate(geom = "text", label = pval, x = 1, y= .45, family = "Avenir", size = 2.8)+
                  labs(x="", y="Dim-1", title = "Periodontitis - Smoker")+
                  MyTheme+
                  theme(aspect.ratio = 1/5,
                        legend.position = "none")
Dim1_PS


# Dim-2 Periodontitis Smokers
pval <- paste0("P=", compare_means(Axis.2 ~ Site, method = "kruskal.test", data = PS)$p.adj)
Dim2_PS <- ggplot(PS, aes(x = Site, y = Axis.2))+
                  geom_boxplot(aes(fill=Site), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                  scale_x_discrete(limits = (levels(HS$Site)))+
                  scale_fill_manual(values = PerioPalette)+
                  ylim(-0.5, 0.5)+
                  annotate(geom = "text", label = pval, x = 2, y = .5, family = "Avenir", size = 2.8, hjust = 0.5)+
                  labs(x="", y="Dim-2", title = "Periodontitis - Smoker")+
                  MyTheme+
                  theme(aspect.ratio = 5/1,
                        legend.position = "none")
Dim2_PS





HNS <- Dim2_HNS + PlotHNS +  plot_spacer() + Dim1_HNS + plot_layout(ncol = 2, nrow = 2, 
                                                      heights = c(5, 1),
                                                      widths = c(1, 5))
HS <- Dim2_HS + PlotHS +  plot_spacer() + Dim1_HS + plot_layout(ncol = 2, nrow = 2, 
                                                                    heights = c(5, 1),
                                                                    widths = c(1, 5))

PNS <- Dim2_PNS + PlotPNS +  plot_spacer() + Dim1_PNS + plot_layout(ncol = 2, nrow = 2, 
                                                                    heights = c(5, 1),
                                                                    widths = c(1, 5))
PS <- Dim2_PS + PlotPS +  plot_spacer() + Dim1_PS + plot_layout(ncol = 2, nrow = 2, 
                                                                heights = c(5, 1),
                                                                widths = c(1, 5))


(HNS | HS)/(PNS|PS)



ggsave(plot = HNS, 
       device = cairo_pdf,
       width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/HNS Figure 3.pdf")

ggsave(plot = HS, 
       device = cairo_pdf,
       width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/HS Figure 3.pdf")


ggsave(plot = PNS, 
       device = cairo_pdf,
       width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/PNS Figure 3.pdf")



ggsave(plot = PS, 
       device = cairo_pdf,
       width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/PS Figure 3.pdf")



