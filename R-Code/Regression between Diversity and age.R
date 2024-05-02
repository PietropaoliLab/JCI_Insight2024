library(phyloseq)
library(ggpubr)
library(ggplot2)
library(patchwork)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Define my subseto of interest
tmp <- subset_samples(psm_perio, Site == "Subgingival" & Smoking == "No")

# Computing Alpha diversity using Observed
alpha <- plot_richness(tmp, x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
Perio <- alpha[["data"]]


# Plotting results
MyTheme <-   theme_bw(12)+
              theme(axis.text = element_text(colour = "black", size = 10),
                    panel.grid = element_blank(),
                    legend.position = "bottom",
                    aspect.ratio = 1/1,
                    legend.key = element_blank(),
                    legend.background=element_blank(),
                    strip.background = element_rect(fill = NA, colour = NA))


MyPalette <- c("#F88B9D", "#8ECEFD")
P <- ggscatter(Perio, x = "AgeYrs", y = "value", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          palette = MyPalette,
          shape = 21,
          alpha = 0.6,
          size = 3,
          fill = "Sex")+
          labs(x="Age (years)", y="Observed ASV diversity", title = "Periodontitis")+
          scale_x_continuous(breaks = seq(30, 75, by = 10))+
          scale_y_continuous(breaks = seq(20, 70, by = 10), limits = c(20, 70))+
          MyTheme+
          stat_cor(aes(color=Sex),method = "pearson")


####################################################################


load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# Define my subseto of interest
rm(tmp)
tmp <- subset_samples(psm_healthy, Site == "Subgingival" & Smoking == "No")

# Computing Alpha diversity using Observed
alpha <- plot_richness(tmp, x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
Healthy <- alpha[["data"]]



H <- ggscatter(Healthy, x = "AgeYrs", y = "value", color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          palette = MyPalette,
          shape = 21,
          size = 3,
          alpha = 0.6,
          fill = "Sex")+
          labs(x="Age (years)", y="Observed ASV diversity", title = "Healthy")+
          scale_x_continuous(breaks = seq(30, 75, by = 10))+
         scale_y_continuous(breaks = seq(20, 70, by = 10), limits = c(20, 70))+
        MyTheme+
          stat_cor(aes(color=Sex),method = "pearson")

H + P

ggsave(plot = H + P, 
       device = cairo_pdf,
       width = 4, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Regression Diversity and age Figure 4.pdf")








