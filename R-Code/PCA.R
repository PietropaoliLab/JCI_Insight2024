library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library("FactoMineR")
library("factoextra")


load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

psm_perio <- tax_glom(psm_perio, taxrank = "Genus")
psm_healthy <- tax_glom(psm_healthy, taxrank = "Genus")


full <- rbind(psmelt(psm_perio), psmelt(psm_healthy))


foo <- full[, c("Genus","Abundance", "Age", "Sex", "Site", "Periodontitis", "Smoking", "BioProject")]
full <- foo[,-1]


full$Age <- as.numeric(full$Age)

# Convert BioProject into numeric
full$BioProject <- as.numeric(as.factor(full$BioProject))


full$Sex <- ifelse(full$Sex=="male", 1, 0)
full$Periodontitis <- ifelse(full$Periodontitis=="Periodontitis", 1, 0)

# 3=NA; 1=Smokers; 2=Non smokers
full$Smoking <- ifelse(is.na(full$Smoking), 3,
                       ifelse(full$Smoking =="Yes", 1, 0))

# Saliva=1; Plaque=2; Subgingivale=3
full$Site <- ifelse(full$Site == "Saliva", 1,
                       ifelse(full$Site =="Plaque", 2, 3))

x <- scale(full, scale = TRUE, center = TRUE)
# PCA
res.pca <- PCA(full,  graph = FALSE)

# Extract eigenvalues/variances
get_eig(res.pca)


# My theme
MyTheme <-  theme_classic(base_size = 10)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

# Control variable colors using their contributions
p0 <- fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping

p0 <- p0 + labs(color="Contribution") +
  MyTheme+
  theme(legend.position = c(.75, .05),
        legend.direction = "horizontal",
        legend.key.height = unit(.2, 'cm'))

# Dimension plot
p <- fviz_screeplot(res.pca, addlabels = TRUE)+
        MyTheme +
        coord_fixed(ylim = c(0, 28), expand = 0)

p[["layers"]][[1]][["aes_params"]][["fill"]] <- "ivory3" # change bar filling
p[["layers"]][[1]][["aes_params"]][["colour"]] <- "black" # No border on barplot
p$layers[[1]]$geom_params$width <- .7 # Change bar width
p[["layers"]][[2]][["aes_params"]][["colour"]] <- "black" # change line color
p[["layers"]][[3]][["aes_params"]][["shape"]] <- 21 # change shape
p[["layers"]][[3]][["aes_params"]][["fill"]] <- "ivory3" # filling gray
p[["layers"]][[3]][["aes_params"]][["colour"]] <- "black" # 
p[["layers"]][[4]][["geom"]][["default_aes"]][["family"]] <- "Avenir"


# Contributions of variables to PC1
p1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)+
        coord_fixed(ylim = c(0, 40), expand = 0)+
        MyTheme+
        theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 0.95))
p1[["labels"]][["title"]] <- "Contribution to Dim-1"
p1[["labels"]][["x"]] <- NULL
p1[["layers"]][[1]][["aes_params"]][["fill"]] <- "ivory3" # change bar filling
p1[["layers"]][[1]][["aes_params"]][["colour"]] <- "black" # No border on barplot

# Contributions of variables to PC2
p2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)+
        coord_fixed(ylim = c(0, 60), expand = 0)+
          MyTheme+
           theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 0.95))
p2[["labels"]][["title"]] <- "Contribution to Dim-2"
p2[["labels"]][["x"]] <- NULL
p2[["layers"]][[1]][["aes_params"]][["fill"]] <- "ivory3" # change bar filling
p2[["layers"]][[1]][["aes_params"]][["colour"]] <- "black" # No border on barplot



library(patchwork)
# PlotPCA <- (p | (p0 / p1 / p2)) + plot_layout(guides = 'collect')

PlotPCA <- (p0 + p) / (p1 + p2)


ggsave(plot = PlotPCA, 
       device = cairo_pdf,
       #width = 8, 
       #height = 11, 
       #units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/PCA.pdf")



PCA_circle <- p0
PCA_dimensions <- p
PCA_Dim1 <- p1
PCA_Dim2 <- p2








