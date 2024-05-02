library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ANCOMBC)
library(dplyr)
library(DESeq2)
library(microbiomeMarker)

set.seed(2022)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

# LEfSe Periodontitis
lefse <- run_lefse(pseq_genus, group = "Sex",
                 taxa_rank = "Genus",
                 kw_cutoff = 1,
                 wilcoxon_cutoff = 1,
                 lda_cutoff = 2)
lefse <- data.frame(lefse@marker_table)
lefse$Genus <- lefse$feature
lefse$Model <- "LEfSe"


# ANCOM-BC Periodontitis
ancom <- ancombc(phyloseq = pseq_genus, 
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

ancom <- data.frame(ancom$res)
colnames(ancom) <- c("beta", "se", "W", "p_val", "padj", "diff_abn")
ancom$Genus <- rownames(ancom)
ancom$Genus <- gsub(".*Genus_", "", ancom$Genus)
ancom$Model = "ANCOM-BC"



# DESeq2 Periodontitis
Deseq <- phyloseq_to_deseq2(pseq_genus, ~ Sex)
# For this dataset need to apply DESeq2 Error in varianceStabilizingTransformation:
 cts <- counts(Deseq)
 geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
 Deseq <- estimateSizeFactors(Deseq, geoMeans=geoMeans)
 Deseq <- DESeq(Deseq, test="Wald", fitType="parametric")
Deseq <- results(Deseq, cooksCutoff = FALSE)
Deseq@listData$Genus <- Deseq@rownames
DeseqRes <- as.data.frame(Deseq@listData)
DeseqRes <- subset(DeseqRes, padj <0.05)
DeseqRes$Genus <- gsub(".*Genus_", "", DeseqRes$Genus)
DeseqRes$Model = "DESeq2"

# Cobine periodontitis results
res <- rbind(DeseqRes[,c("Genus", "padj", "Model")], ancom[,c("Genus", "padj", "Model")])
res <- rbind(res, lefse[,c("Genus", "padj", "Model")])
res <- subset(res, padj <0.05)
res$Group = "Periodontitis"
table(res$Model)
table(res$Genus, res$Model)

Perio <- res
rm(cts, ancom, Deseq, DeseqRes, geoMeans, lefse, pseq_genus)





# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")

# LEfSe Healthy
lefse <- run_lefse(pseq_genus, group = "Sex",
                   taxa_rank = "Genus",
                   kw_cutoff = 1,
                   wilcoxon_cutof = 1,
                   lda_cutoff = 2)
lefse <- data.frame(lefse@marker_table)
lefse$Genus <- lefse$feature
lefse$Model <- "LEfSe"

# ANCOM-BC Healthy
ancom <- ancombc(phyloseq = pseq_genus, 
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

ancom <- data.frame(ancom$res)
colnames(ancom) <- c("beta", "se", "W", "p_val", "padj", "diff_abn")
ancom$Genus <- rownames(ancom)
ancom$Genus <- gsub(".*Genus_", "", ancom$Genus)
ancom$Model = "ANCOM-BC"



# DESeq2 Healthy
Deseq <- phyloseq_to_deseq2(pseq_genus, ~ Sex)
# For this dataset need to apply DESeq2 Error in varianceStabilizingTransformation:
cts <- counts(Deseq)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
Deseq <- estimateSizeFactors(Deseq, geoMeans=geoMeans)
Deseq <- DESeq(Deseq, test="Wald", fitType="parametric")
Deseq <- results(Deseq, cooksCutoff = FALSE)
Deseq@listData$Genus <- Deseq@rownames
DeseqRes <- as.data.frame(Deseq@listData)
DeseqRes$Genus <- gsub(".*Genus_", "", DeseqRes$Genus)
DeseqRes$Model = "DESeq2"
res <- rbind(DeseqRes[,c("Genus", "padj", "Model")], ancom[,c("Genus", "padj", "Model")])
res <- rbind(res, lefse[,c("Genus", "padj", "Model")])
res <- subset(res, padj <0.05)
res$Group = "Healthy"
Healthy <- res

table(res$Model)
table(res$Genus, res$Model)

res <- rbind(Healthy, Perio)


rm(cts, ancom, Deseq, DeseqRes, geoMeans, lefse, pseq_genus, psm_healthy, psm_perio)





table(Perio$Model)
table(Perio$Genus, Perio$Model)





table(Healthy$Model)
table(Healthy$Genus, Healthy$Model)

