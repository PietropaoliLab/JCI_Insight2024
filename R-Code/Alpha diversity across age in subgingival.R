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
        legend.position = "none",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

set.seed(101)
p1 <- plot_richness(tmp, x="AgeCat", measures="Shannon", color = "Sex")

# My Shannon Index
shannon <- p1$data

p <- ggline(subset(shannon, Site == "Subgingival"), x = "AgeCat", y = "value", 
       add = c("mean_se"),
       facet.by = c("Periodontitis"),
       color = "Sex", palette = c("#f2acda", "#6c9bd9"))+
  stat_compare_means(label.y = 1.7)+
  labs(x="Age", y="Subgingival diversity (Shannon)")+
  MyTheme

# Dot size
p[["layers"]][[3]][["aes_params"]][["size"]] <- 2.6

# Line size
p[["layers"]][[2]][["aes_params"]][["size"]] <- .75

# Text of Statistic informations
p[["layers"]][[4]][["geom"]][["default_aes"]][["family"]] <- "Avenir"
p[["layers"]][[4]][["geom"]][["default_aes"]][["size"]] <- 3.5
p



# stat_compare_means(aes(group = Sex), label = "p.format")  
  # Not significance in comparisons F/M

