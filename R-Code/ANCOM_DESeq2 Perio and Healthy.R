library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ANCOMBC)
library(dplyr)
library(DESeq2)

# ANCOM
# The ANCOM procedure compares the relative abundance of a taxon between two ecosystems by computing Aitchison’s [5] 
# log-ratio of abundance of each taxon relative to the abundance of all remaining taxa one at a time. Thus, if there
# are “m” taxa, then for each taxon it performs “m-1” tests and the significance of each test is determined using 
# the Benjamini-Hochberg procedure that controls for FDR at 0.05. For each taxon, ANCOM counts the number of tests 
# among the m-1 tests that are rejected. Thus for each taxon, ANCOM obtains a count random variable W that 
# represents the number of nulls among the m-1 tests that are rejected. ANCOM determines the final significance of 
# a taxon by using the empirical distribution of W. To deal with zero counts, we use an arbitrary pseudo count 
# value of 0.001. For a more detailed description of ANCOM, we refer the reader to Mandal et al. 

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")

Plaque <- subset_samples(pseq_genus, Site == "Plaque")
Saliva <- subset_samples(pseq_genus, Site == "Saliva")
Subgingival <- subset_samples(pseq_genus, Site == "Subgingival")



# ####################################

PlaqueOut = ancombc(phyloseq = Plaque, 
                    formula = "Sex", 
                    p_adj_method = "fdr", 
                    zero_cut = 0.90, # by default prevalence filter of 10% is applied
                    lib_cut = 0, 
                    group = "Sex", 
                    struc_zero = TRUE, 
                    neg_lb = TRUE, 
                    tol = 1e-5, 
                    max_iter = 100, 
                    conserve = TRUE, 
                    alpha = 0.05, 
                    global = TRUE)

SalivaOut = ancombc(phyloseq = Saliva, 
                    formula = "Sex", 
                    p_adj_method = "fdr", 
                    zero_cut = 0.90, # by default prevalence filter of 10% is applied
                    lib_cut = 0, 
                    group = "Sex", 
                    struc_zero = TRUE, 
                    neg_lb = TRUE, 
                    tol = 1e-5, 
                    max_iter = 100, 
                    conserve = TRUE, 
                    alpha = 0.05, 
                    global = TRUE)

SubgingivalOut = ancombc(phyloseq = Subgingival, 
                    formula = "Sex", 
                    p_adj_method = "fdr", 
                    zero_cut = 0.90, # by default prevalence filter of 10% is applied
                    lib_cut = 0, 
                    group = "Sex", 
                    struc_zero = TRUE, 
                    neg_lb = TRUE, 
                    tol = 1e-5, 
                    max_iter = 100, 
                    conserve = TRUE, 
                    alpha = 0.05, 
                    global = TRUE)

PlaqueRes = PlaqueOut$res
SalivaRes = SalivaOut$res
SubgingivalRes = SubgingivalOut$res


PerioRes <- rbind(data.frame(PlaqueRes$diff_abn, "Site" = "Plaque", "p.adj" = PlaqueRes$q_val),
             data.frame(SalivaRes$diff_abn, "Site" = "Saliva", "p.adj" = SalivaRes$q_val))
PerioRes <- rbind(PerioRes, data.frame(SubgingivalRes$diff_abn, "Site" = "Subgingival", "p.adj" = SubgingivalRes$q_val))

rm(PlaqueRes, SalivaRes, SubgingivalRes, Plaque, Saliva, Subgingival, PlaqueOut, SalivaOut, SubgingivalOut, psm_perio, pseq_genus)


##########################################################################

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)

Plaque <- subset_samples(pseq_genus, Site == "Plaque")
Saliva <- subset_samples(pseq_genus, Site == "Saliva")
Subgingival <- subset_samples(pseq_genus, Site == "Subgingival")



# ####################################

PlaqueOut = ancombc(phyloseq = Plaque, 
                    formula = "Sex", 
                    p_adj_method = "fdr", 
                    zero_cut = 0.90, # by default prevalence filter of 10% is applied
                    lib_cut = 0, 
                    group = "Sex", 
                    struc_zero = TRUE, 
                    neg_lb = TRUE, 
                    tol = 1e-5, 
                    max_iter = 100, 
                    conserve = TRUE, 
                    alpha = 0.05, 
                    global = TRUE)

SalivaOut = ancombc(phyloseq = Saliva, 
                    formula = "Sex", 
                    p_adj_method = "fdr", 
                    zero_cut = 0.90, # by default prevalence filter of 10% is applied
                    lib_cut = 0, 
                    group = "Sex", 
                    struc_zero = TRUE, 
                    neg_lb = TRUE, 
                    tol = 1e-5, 
                    max_iter = 100, 
                    conserve = TRUE, 
                    alpha = 0.05, 
                    global = TRUE)

SubgingivalOut = ancombc(phyloseq = Subgingival, 
                         formula = "Sex", 
                         p_adj_method = "fdr", 
                         zero_cut = 0.90, # by default prevalence filter of 10% is applied
                         lib_cut = 0, 
                         group = "Sex", 
                         struc_zero = TRUE, 
                         neg_lb = TRUE, 
                         tol = 1e-5, 
                         max_iter = 100, 
                         conserve = TRUE, 
                         alpha = 0.05, 
                         global = TRUE)

PlaqueRes = PlaqueOut$res
SalivaRes = SalivaOut$res
SubgingivalRes = SubgingivalOut$res

# ####################################


HealthyRes <- rbind(data.frame(PlaqueRes$diff_abn, "Site" = "Plaque", "p.adj" = PlaqueRes$q_val),
                  data.frame(SalivaRes$diff_abn, "Site" = "Saliva", "p.adj" = SalivaRes$q_val))
HealthyRes <- rbind(HealthyRes, data.frame(SubgingivalRes$diff_abn, "Site" = "Subgingival", "p.adj" = SubgingivalRes$q_val))



rm(PlaqueRes, SalivaRes, SubgingivalRes, Plaque, Saliva, Subgingival, PlaqueOut, SalivaOut, SubgingivalOut, psm_healthy, pseq_genus)



ANCOMData <- rbind(data.frame(HealthyRes, "Group" = "Healthy"), 
                   data.frame(PerioRes, "Group" = "Periodontitis"))
ANCOMData$Genus <- rownames(ANCOMData)

rm(HealthyRes, PerioRes)

####################################################################################################
# Computing
# DESEq2 

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")

# Subsetting my data
Plaque <- subset_samples(pseq_genus, Site == "Plaque")
Saliva <- subset_samples(pseq_genus, Site == "Saliva")
Subgingival <- subset_samples(pseq_genus, Site == "Subgingival")

# Set my cutoff
alpha <- 0.01

# Deseq2 operations

# Plaque
Plaquediagdds <- phyloseq_to_deseq2(Plaque, ~ Sex)
Plaquediagdds <- DESeq(Plaquediagdds, test="Wald", fitType="parametric")
PlaqueRes <- results(Plaquediagdds, cooksCutoff = FALSE)
PlaqueRes@listData$Genus <- PlaqueRes@rownames
PlaqueRes@listData$Site = "Plaque"
PlaqueRes@listData$Group = "Periodontitis"


# Saliva
Salivadiagdds <- phyloseq_to_deseq2(Saliva, ~ Sex)
# For this dataset need to apply DESeq2 Error in varianceStabilizingTransformation:
cts <- counts(Salivadiagdds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
Salivadiagdds <- estimateSizeFactors(Salivadiagdds, geoMeans=geoMeans)
Salivadiagdds <- DESeq(Salivadiagdds, test="Wald", fitType="parametric")
SalivaRes <- results(Salivadiagdds, cooksCutoff = FALSE)
SalivaRes@listData$Genus <- SalivaRes@rownames
SalivaRes@listData$Site = "Saliva"
SalivaRes@listData$Group = "Periodontitis"


# Subgingival
Subgingivaldiagdds <- phyloseq_to_deseq2(Subgingival, ~ Sex)
Subgingivaldiagdds <- DESeq(Subgingivaldiagdds, test="Wald", fitType="parametric")
SubgingivalRes <- results(Subgingivaldiagdds, cooksCutoff = FALSE)
SubgingivalRes@listData$Genus <- SubgingivalRes@rownames
SubgingivalRes@listData$Site = "Subgingival"
SubgingivalRes@listData$Group = "Periodontitis"

DESeqPerio <- rbind(as.data.frame(PlaqueRes), as.data.frame(SalivaRes))
DESeqPerio <- rbind(DESeqPerio, as.data.frame(SubgingivalRes))

rm("alpha", "cts",  "diagdds", "geoMeans", 
"Plaque", "Plaquediagdds", "PlaqueRes", "pseq_genus", "psm_perio", 
"Saliva", "Salivadiagdds", "SalivaRes", "Subgingival", "Subgingivaldiagdds", 
"SubgingivalRes")


######
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")

# Subsetting my data
Plaque <- subset_samples(pseq_genus, Site == "Plaque")
Saliva <- subset_samples(pseq_genus, Site == "Saliva")
Subgingival <- subset_samples(pseq_genus, Site == "Subgingival")

# Set my cutoff
alpha <- 0.01

# Deseq2 operations

# Plaque
Plaquediagdds <- phyloseq_to_deseq2(Plaque, ~ Sex)
Plaquediagdds <- DESeq(Plaquediagdds, test="Wald", fitType="parametric")
PlaqueRes <- results(Plaquediagdds, cooksCutoff = FALSE)
PlaqueRes@listData$Genus <- PlaqueRes@rownames
PlaqueRes@listData$Site = "Plaque"
PlaqueRes@listData$Group = "Healthy"


# Saliva
Salivadiagdds <- phyloseq_to_deseq2(Saliva, ~ Sex)
# For this dataset need to apply DESeq2 Error in varianceStabilizingTransformation:
# cts <- counts(Salivadiagdds)
# geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
# Salivadiagdds <- estimateSizeFactors(Salivadiagdds, geoMeans=geoMeans)
Salivadiagdds <- DESeq(Salivadiagdds, test="Wald", fitType="parametric")
SalivaRes <- results(Salivadiagdds, cooksCutoff = FALSE)
SalivaRes@listData$Genus <- SalivaRes@rownames
SalivaRes@listData$Site = "Saliva"
SalivaRes@listData$Group = "Healthy"


# Subgingival
Subgingivaldiagdds <- phyloseq_to_deseq2(Subgingival, ~ Sex)
Subgingivaldiagdds <- DESeq(Subgingivaldiagdds, test="Wald", fitType="parametric")
SubgingivalRes <- results(Subgingivaldiagdds, cooksCutoff = FALSE)
SubgingivalRes@listData$Genus <- SubgingivalRes@rownames
SubgingivalRes@listData$Site = "Subgingival"
SubgingivalRes@listData$Group = "Healthy"

DESeqHealthy <- rbind(as.data.frame(PlaqueRes), as.data.frame(SalivaRes))
DESeqHealthy <- rbind(DESeqHealthy, as.data.frame(SubgingivalRes))

rm("alpha", "cts",  "diagdds", "geoMeans", 
   "Plaque", "Plaquediagdds", "PlaqueRes", "pseq_genus", "psm_healthy", 
   "Saliva", "Salivadiagdds", "SalivaRes", "Subgingival", "Subgingivaldiagdds", 
   "SubgingivalRes")


DESeqData <- rbind(DESeqPerio, DESeqHealthy)
rm(DESeqPerio, DESeqHealthy)


gc()


colnames(ANCOMData) <- c("Selected", "Site", "padj", "Group", "Genus")

foo <- full_join(ANCOMData[,c("Site", "padj", "Group", "Genus")], DESeqData[,c("Site", "padj", "Group", "Genus")], by = "Genus",
                 suffix = c("_ANCOM", "_DESeq"))


PlaqueHealthy <- subset(foo, Site_ANCOM == "Plaque" & Site_DESeq == "Plaque" & Group_ANCOM == "Healthy" & Group_DESeq == "Healthy")
PlaqueHealthy$Group <- "Healthy"
PlaqueHealthy$Site <- "Plaque"
PlaqueHealthy <- PlaqueHealthy[,c(2,4, 6, 8, 9)]

PlaquePerio <- subset(foo, Site_ANCOM == "Plaque" & Site_DESeq == "Plaque" & Group_ANCOM == "Periodontitis" & Group_DESeq == "Periodontitis")
PlaquePerio$Group <- "Periodontitis"
PlaquePerio$Site <- "Plaque"
PlaquePerio <- PlaquePerio[,c(2,4, 6, 8, 9)]


SalivaHealthy <- subset(foo, Site_ANCOM == "Saliva" & Site_DESeq == "Saliva" & Group_ANCOM == "Healthy" & Group_DESeq == "Healthy")
SalivaHealthy$Group <- "Healthy"
SalivaHealthy$Site <- "Saliva"
SalivaHealthy <- SalivaHealthy[,c(2,4, 6, 8, 9)]

SalivaPerio <- subset(foo, Site_ANCOM == "Saliva" & Site_DESeq == "Saliva" & Group_ANCOM == "Periodontitis" & Group_DESeq == "Periodontitis")
SalivaPerio$Group <- "Periodontitis"
SalivaPerio$Site <- "Saliva"
SalivaPerio <- SalivaPerio[,c(2,4, 6, 8, 9)]


SubgingivalHealthy <- subset(foo, Site_ANCOM == "Subgingival" & Site_DESeq == "Subgingival" & Group_ANCOM == "Healthy" & Group_DESeq == "Healthy")
SubgingivalHealthy$Group <- "Healthy"
SubgingivalHealthy$Site <- "Subgingival"
SubgingivalHealthy <- SubgingivalHealthy[,c(2,4, 6, 8, 9)]

SubgingivalPerio <- subset(foo, Site_ANCOM == "Subgingival" & Site_DESeq == "Subgingival" & Group_ANCOM == "Periodontitis" & Group_DESeq == "Periodontitis")
SubgingivalPerio$Group <- "Periodontitis"
SubgingivalPerio$Site <- "Subgingival"
SubgingivalPerio <- SubgingivalPerio[,c(2,4, 6, 8, 9)]


library(dplyr)

DataForPlot <- bind_rows(PlaqueHealthy, PlaquePerio, SalivaHealthy, SalivaPerio, SubgingivalHealthy, SubgingivalPerio)
DataForPlot <- DataForPlot[,c(2, 1, 3:5)]
DataForPlot$Select <- ifelse(DataForPlot$padj_ANCOM <0.05 & DataForPlot$padj_DESeq <0.05, "Both FDR<0.05", NA)
DataForPlot$Genus <- gsub(".*Genus_", "", DataForPlot$Genus)
DataForPlot$MyLabs <- ifelse(DataForPlot$Select == "Both FDR<0.05", DataForPlot$Genus, NA)



rm(PlaqueHealthy, PlaquePerio, SalivaHealthy, SalivaPerio, SubgingivalHealthy, SubgingivalPerio)

library(ggplot2)
library(ggrepel)
DataForPlot <- DataForPlot[order(DataForPlot$padj_DESeq, decreasing=TRUE), ]

MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank())

ggplot(DataForPlot, aes(x=padj_ANCOM, y=padj_DESeq, group = Select))+
  geom_point(aes(color = Select, shape = Site), size = 3.5)+
  scale_color_manual(values = c("#ca0020"),  na.value = "#d9d9d9") +
  #scale_shape_manual(values=c(15, 16 , 17))+
  #coord_fixed(xlim = c(0, .45), ylim = c(0, .45))+
  geom_label_repel(aes(label = MyLabs), size = 3, fill="white", seed = 101,
                   box.padding = unit(1.5, "lines"),
                   point.padding = unit(0.3, "lines"))+
  labs(x = "P value ANCOM-BC", y= "P value DESeq2")+
  facet_wrap(~Group)+
  scale_y_reverse() +
  scale_x_reverse() +
  MyTheme+
  theme(legend.position = "right",
        legend.direction = "vertical")











