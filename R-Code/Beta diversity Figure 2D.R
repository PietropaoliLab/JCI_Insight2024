library("microbial")
library(phyloseq)
library(microbiome)
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
        legend.position = "right",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

# Healthy
Beta_Healthy <- plotbeta(psm_healthy, group="Site", 
                distance = "bray", 
                ellipse = FALSE, 
                color = c("#ccece6", "#1c9099", "#006d2c"))+
                    scale_x_continuous(breaks = c(seq(-0.4, 0.4, 0.2)))+
                    MyTheme + 
                    coord_fixed(xlim = c(-0.43, 0.43)) + 
                    labs(title = "Healthy", color = "Site")


# add alpha               
Beta_Healthy$layers[[1]]$aes_params$alpha <- .8
Beta_Healthy$layers[[1]]$aes_params$size <- 3
Beta_Healthy

# Periodontitis
Beta_Perio <- plotbeta(psm_perio, group="Site", 
                         distance = "bray", 
                         ellipse = FALSE, 
                         color = c("#fc9272", "#531C46", "#B9263A"))+
                         scale_x_continuous(breaks = c(seq(-0.4, 0.4, 0.2)))+
                            MyTheme + 
                            coord_fixed(xlim = c(-0.50, 0.43)) + 
                            labs(title = "Periodontitis", color = "Site")
# add alpha               
Beta_Perio$layers[[1]]$aes_params$alpha <- .8
Beta_Perio$layers[[1]]$aes_params$size <- 3


Beta_Perio

# Combining plots
Beta <- Beta_Healthy + Beta_Perio + plot_layout(ncol = 1)

Beta

ggsave(plot = Beta, 
       device = cairo_pdf,
       width = 4, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Beta_diversity_Figure 2D.pdf")


# ===================================== #

#     Compute P-values for Brays
#         Kruscal wallis
# ===================================== #

DatasetHealthy <- Beta_Healthy$data
DatasetPerio <- Beta_Perio$data

# Healthy
compare_means(Axis.1 ~ Site, method = "kruskal.test", data = DatasetHealthy)

# Periodontitis
compare_means(Axis.1 ~ Site, method = "kruskal.test", data = DatasetPerio)



HealthyPalette <- c("#ccece6", "#1c9099", "#006d2c")
Dim1_healthy <- ggplot(DatasetHealthy, aes(x = group, y = Axis.1))+
                      geom_boxplot(aes(fill=group), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                      scale_x_discrete(limits = rev(levels(DatasetHealthy$group)))+
                      scale_fill_manual(values = HealthyPalette)+
                      scale_y_continuous(breaks = c(seq(-0.4, 0.4, 0.2)))+
                      coord_fixed(ylim = c(-0.43, 0.43))+
                            coord_flip()+
                            labs(x="", y="Dim-1")+
                            MyTheme+
                            theme(aspect.ratio = 1/5,
                                  legend.position = "none")

Dim2_healthy <- ggplot(DatasetHealthy, aes(x = group, y = Axis.2))+
                    geom_boxplot(aes(fill=group), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                    scale_fill_manual(values = HealthyPalette)+
                    labs(x="", y="Dim-2")+
                    MyTheme+
                    theme(aspect.ratio = 5/1,
                          legend.position = "none")
compare_means(Axis.1 ~ Site, method = "kruskal.test", data = DatasetHealthy)
compare_means(Axis.2 ~ Site, method = "kruskal.test", data = DatasetHealthy)




PerioPalette <- c("#fc9272", "#531C46", "#B9263A")
Dim1_perio <- ggplot(DatasetPerio, aes(x = group, y = Axis.1))+
                        geom_boxplot(aes(fill=group), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                        scale_x_discrete(limits = rev(levels(DatasetPerio$group)))+
                        scale_fill_manual(values = PerioPalette)+
                        scale_y_continuous(breaks = c(seq(-0.4, 0.4, 0.2)))+
                        coord_fixed(ylim = c(-0.43, 0.43))+
                        coord_flip()+
                        labs(x="", y="Dim-1")+
                        MyTheme+
                        theme(aspect.ratio = 1/5,
                              legend.position = "none")

Dim2_perio <- ggplot(DatasetPerio, aes(x = group, y = Axis.2))+
                      geom_boxplot(aes(fill=group), width = .75, color = "#ceccc1", outlier.shape = NA, alpha = .8)+
                      scale_fill_manual(values = PerioPalette)+
                      labs(x="", y="Dim-2")+
                      MyTheme+
                      theme(aspect.ratio = 5/1,
                            legend.position = "none")

compare_means(Axis.1 ~ Site, method = "kruskal.test", data = DatasetPerio)
compare_means(Axis.2 ~ Site, method = "kruskal.test", data = DatasetPerio)

ggsave(plot = Dim1_perio, 
       device = cairo_pdf,
       #width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Dim1_perio.pdf")
ggsave(plot = Dim2_perio, 
       device = cairo_pdf,
       #width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Dim2_perio.pdf")


ggsave(plot = Dim1_healthy, 
       device = cairo_pdf,
       #width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Dim1_healthy.pdf")
ggsave(plot = Dim2_healthy, 
       device = cairo_pdf,
       #width = 5, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Dim2_healthy.pdf")



