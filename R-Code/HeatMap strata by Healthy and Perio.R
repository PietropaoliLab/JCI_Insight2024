##########################################################################
# Generating HeatMap with Predictions
library(ComplexHeatmap)
library(reshape2)
library(microbiome)
library(dplyr)
library(circlize)



# Loading datasets
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
healthy <- phyloseq::tax_glom(psm_healthy, taxrank = "Phylum")
perio <- phyloseq::tax_glom(psm_perio, taxrank = "Phylum")

# Relative abundance
healthy <- psm_healthy %>%
                          aggregate_taxa(level = "Phylum") %>%  
                          microbiome::transform(transform = "compositional")
perio <- psm_perio %>%
                        aggregate_taxa(level = "Phylum") %>%  
                        microbiome::transform(transform = "compositional")


# Melting Phyloseq Objects
healthy <- psmelt(healthy)
perio <- psmelt(perio)

# Sorting
healthy <- with(healthy, healthy[order(Site, Sex),])
perio <- with(perio, perio[order(Site, Sex),])

healthy$Sex <- factor(healthy$Sex, levels = c("female", "male"), labels = c("F", "M"))
perio$Sex <- factor(perio$Sex, levels = c("female", "male"), labels = c("F", "M"))



HealthyMatrix <- dcast(healthy,  Phylum ~ Run, value.var = "Abundance", fun.aggregate=sum)
rownames(HealthyMatrix) <- HealthyMatrix$Phylum
HealthyMatrix <- HealthyMatrix[-1]

# Remove Phyla with relative abundance = 0%
ToRemove <- c("Acidobacteriota", "Euryarchaeota")
HealthyMatrix <- HealthyMatrix[!rownames(HealthyMatrix) %in% ToRemove, ]  # ! is logical negation
# Scaling by rows
HealthyMatrix <- t(scale(t(HealthyMatrix)))


# Define annotations for the heatmaps
Annotations <- healthy[,c("Run", "Site", "Sex")]
Annotations <- unique(Annotations)

# defining my palette
col_fun = colorRamp2(c(0, 0.3, 0.6), c("#f7f7f7", "#fddbc7", "#b2182b"))
col_fun(seq(-3, 3))


# Save PDF
pdf(file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Heatmap Healthy.pdf",
    width = 20, height = 3.5)
# Drawing Heatmap

HPhealthy <- Heatmap(HealthyMatrix, name = "Relative abundance\nz-score",
                        cluster_columns = FALSE,
                        col = col_fun,
                        rect_gp = gpar(col = "gray80", lwd = 0.5),
                        row_names_gp = gpar(fontsize = 10),
                        column_split = Annotations$Site,
                        top_annotation = HeatmapAnnotation(Sex = Annotations$Sex,  which = 'col', 
                                           col = list('Sex' = c('F' = 'pink', 'M' = 'royalblue'))),
                        column_gap = unit(c(3, 3, 3), "mm"),
                        #border = TRUE,
                        # bottom_annotation = HeatmapAnnotation(Sex = SalivaAnnotation$Sex), 
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        show_row_names = TRUE,
                        row_names_side = "left",
                        show_column_names = FALSE)
draw(HPhealthy)
dev.off()
#####################################################################################################à

PerioMatrix <- dcast(perio,  Phylum ~ Run, value.var = "Abundance", fun.aggregate=sum)
rownames(PerioMatrix) <- PerioMatrix$Phylum
PerioMatrix <- PerioMatrix[-1]

# Remove Phyla with relative abundance = 0%
# ToRemove <- c("Chloroflexi", "Acidobacteriota", "Euryarchaeota", "Deinococcota", "Verrucomicrobiota")
# HealthyMatrix <- HealthyMatrix[!rownames(HealthyMatrix) %in% ToRemove, ]  # ! is logical negation

# Scaling by rows
PerioMatrix <- t(scale(t(PerioMatrix)))

# Define annotations for the heatmaps
Annotations <- perio[,c("Run", "Site", "Sex")]
Annotations <- unique(Annotations)

# defining my palette
col_fun = colorRamp2(c(0, 0.3, 0.6), c("#f7f7f7", "#fddbc7", "#b2182b"))
col_fun(seq(-3, 3))

# Save PDF
pdf(file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Heatmap Perio.pdf",
    width = 20, height = 3.5)
# Drawing Heatmap
HMperio <- Heatmap(PerioMatrix, name = "Relative abundance\nz-score",
                        cluster_columns = FALSE,
                        col = col_fun,
                        rect_gp = gpar(col = "gray80", lwd = 0.2),
                        row_names_gp = gpar(fontsize = 10),
                        column_split = Annotations$Site,
                        bottom_annotation = HeatmapAnnotation(Sex = Annotations$Sex,  which = 'col', 
                                           col = list('Sex' = c('F' = 'pink', 'M' = 'royalblue'))),
                        column_gap = unit(c(3, 3, 3), "mm"),
                        #border = TRUE,
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        show_row_names = TRUE,
                        row_names_side = "left",
                        show_column_names = FALSE)

draw(HMperio)
dev.off()



