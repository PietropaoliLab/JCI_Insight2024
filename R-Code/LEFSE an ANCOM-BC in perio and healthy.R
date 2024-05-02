library(phyloseq)
library(microbiomeMarker)
library(ANCOMBC)
library(dplyr)



# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

ps <- tax_glom(psm_perio, taxrank = "Genus")
MySite <- c("Saliva", "Plaque", "Subgingival")
res <- data.frame()

for (i in MySite) {
  foo <- subset_samples(ps, Site == i)
  
  tmp <- run_lefse(foo, group = "Sex",
                 taxa_rank = "Genus",
                 kw_cutoff = 1,
                 wilcoxon_cutoff	= 1)
  tmp <- data.frame(tmp@marker_table)
  tmp$Site = i
  tmp$Group = "Periodontitis"
  res <- rbind(res, tmp)
  rm(foo, tmp)
  
}
rm(ps)



# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

ps <- tax_glom(psm_healthy, taxrank = "Genus")

for (i in MySite) {
  foo <- subset_samples(ps, Site == i)
  
  tmp <- run_lefse(foo, group = "Sex",
                   taxa_rank = "Genus",
                   kw_cutoff = 1,
                   wilcoxon_cutoff	= 1)
  tmp <- data.frame(tmp@marker_table)
  tmp$Site = i
  tmp$Group = "Healthy"
  res <- rbind(res, tmp)
  rm(foo, tmp)
  
}
LefseData <- res
rm(ps, MySite, i, res)

###################################################################################

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

colnames(ANCOMData) <- c("Selected", "Site", "padj", "Group", "Genus")

####################################################################################################

# Combining dataframes
ANCOMData$Genus <- gsub(".*Genus_", "", ANCOMData$Genus)

colnames(LefseData) <- c("Genus", "Sex","ef_lda", "pvalue", "padj", "Site", "Group")

foo <- full_join(ANCOMData[,c("Site", "padj", "Group", "Genus")], LefseData[,c("Site", "padj", "Group", "Genus")], by = "Genus",
                 suffix = c("_ANCOM", "_Lefse"))


PlaqueHealthy <- subset(foo, Site_ANCOM == "Plaque" & Site_Lefse == "Plaque" & Group_ANCOM == "Healthy" & Group_Lefse == "Healthy")
PlaqueHealthy$Group <- "Healthy"
PlaqueHealthy$Site <- "Plaque"
PlaqueHealthy <- PlaqueHealthy[,c(2,4, 6, 8, 9)]

PlaquePerio <- subset(foo, Site_ANCOM == "Plaque" & Site_Lefse == "Plaque" & Group_ANCOM == "Periodontitis" & Group_Lefse == "Periodontitis")
PlaquePerio$Group <- "Periodontitis"
PlaquePerio$Site <- "Plaque"
PlaquePerio <- PlaquePerio[,c(2,4, 6, 8, 9)]

SalivaHealthy <- subset(foo, Site_ANCOM == "Saliva" & Site_Lefse == "Saliva" & Group_ANCOM == "Healthy" & Group_Lefse == "Healthy")
SalivaHealthy$Group <- "Healthy"
SalivaHealthy$Site <- "Saliva"
SalivaHealthy <- SalivaHealthy[,c(2,4, 6, 8, 9)]

SalivaPerio <- subset(foo, Site_ANCOM == "Saliva" & Site_Lefse == "Saliva" & Group_ANCOM == "Periodontitis" & Group_Lefse == "Periodontitis")
SalivaPerio$Group <- "Periodontitis"
SalivaPerio$Site <- "Saliva"
SalivaPerio <- SalivaPerio[,c(2,4, 6, 8, 9)]

SubgingivalHealthy <- subset(foo, Site_ANCOM == "Subgingival" & Site_Lefse == "Subgingival" & Group_ANCOM == "Healthy" & Group_Lefse == "Healthy")
SubgingivalHealthy$Group <- "Healthy"
SubgingivalHealthy$Site <- "Subgingival"
SubgingivalHealthy <- SubgingivalHealthy[,c(2,4, 6, 8, 9)]

SubgingivalPerio <- subset(foo, Site_ANCOM == "Subgingival" & Site_Lefse == "Subgingival" & Group_ANCOM == "Periodontitis" & Group_Lefse == "Periodontitis")
SubgingivalPerio$Group <- "Periodontitis"
SubgingivalPerio$Site <- "Subgingival"
SubgingivalPerio <- SubgingivalPerio[,c(2,4, 6, 8, 9)]



# Combining dataframes and prepare for ggplot
DataForPlot <- bind_rows(PlaqueHealthy, PlaquePerio, SalivaHealthy, SalivaPerio, SubgingivalHealthy, SubgingivalPerio)
DataForPlot <- DataForPlot[,c(2, 1, 3:5)]
DataForPlot$Select <- ifelse(DataForPlot$padj_ANCOM <0.05 & DataForPlot$padj_Lefse <0.05, "Both FDR<0.05", NA)
DataForPlot$MyLabs <- ifelse(DataForPlot$Select == "Both FDR<0.05", DataForPlot$Genus, NA)

# some cleening
rm(PlaqueHealthy, PlaquePerio, SalivaHealthy, SalivaPerio, SubgingivalHealthy, SubgingivalPerio, foo)

# plotting
library(ggplot2)
library(ggrepel)
DataForPlot <- DataForPlot[order(DataForPlot$padj_Lefse, decreasing=TRUE), ]

MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank())

ggplot(DataForPlot, aes(x=padj_ANCOM, y=padj_Lefse, group = Select))+
  geom_point(aes(color = Select, shape = Site), size = 3.5)+
  scale_color_manual(values = c("#ca0020"),  na.value = "#d9d9d9") +
  #scale_shape_manual(values=c(15, 16 , 17))+
  #coord_fixed(xlim = c(0, .45), ylim = c(0, .45))+
  geom_label_repel(aes(label = MyLabs), size = 2.5, fill="white", seed = 101,
                   box.padding = 0.3,
                   point.padding = 0.3,
                   max.overlaps = 15)+
  labs(x = "P value ANCOM-BC", y= "P value LEfSe")+
  facet_wrap(~Group)+
  scale_y_reverse() +
  scale_x_reverse() +
  MyTheme+
  theme(legend.position = "right",
        legend.direction = "vertical")







