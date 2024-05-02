library(phyloseq)
library(dplyr)
library(microbiomeMarker)
library(ggpubr)
library(ggplot2)

# Define my storing list
MLresults <- list()

# My palette
MyPalette <- c("#F88B9D", "#8ECEFD")

# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))


load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# my results df
results <- data.frame()


# Define my subseto of interest
tmp <- subset_samples(psm_perio, Site == "Subgingival" & Smoking == "No")

# Run t-test for DA
foo <- run_simple_stat(tmp,
                       taxa_rank = "Genus",
                       transform = "log10p",
                       norm = "TSS", #  total sum scaling, also referred to as "relative abundance", the abundances were normalized by dividing the corresponding sample library size.
                       method = "t.test",
                       # p_adjust = "fdr",
                       group = "Sex")

results <- data.frame(marker_table(foo))
results$Model = "t-test"
results$Smoking = "Non-smoking"

# Apply FDR
results$padj <- p.adjust(results$pvalue, method = "fdr")


# Plotting results
PlotPerioEnrichedGenera <- ggplot(results, aes(x=reorder(feature, ef_diff_mean), y = ef_diff_mean, fill = enrich_group))+
                                    geom_hline(yintercept = 0, size = .1)+
                                    geom_col(width = .75)+
                                    scale_fill_manual(values = MyPalette)+
                                    labs(x="", y="Effect size estimated mean")+
                                    coord_flip()+
                                      MyTheme
PlotPerioEnrichedGenera

# Store my enriched genera plot 
MLresults[["PlotPerioEnrichedGenera"]] <- PlotPerioEnrichedGenera

####################################
#
##  Microgenderome Index - MGI
#
####################################
tmp <- tax_glom(tmp, taxrank = "Genus")
df_alpha <- tmp
tmp <- microbiome::transform(tmp, transform = "log10p")
tmp <- microbiome::transform(tmp, transform = "compositional")


# Subset Genera of interest
MaleGenus <- na.omit(ifelse(results$ef_diff_mean < 0, results$feature, NA))
FemaleGenus <- na.omit(ifelse(results$ef_diff_mean > 0, results$feature, NA))

# Prepare data for MGI
df <- subset_taxa(tmp, Genus %in% results$feature)
df <- data.frame(otu_table(df))
colnames(df) <- gsub(".*Genus_", "", colnames(df))

# Sums relative abundance of enriched genus in male or in female
df$MaleEnriched <- rowSums(df[MaleGenus])
df$FemaleEnriched <- rowSums(df[FemaleGenus])

# Compute MGI
df$MGI <- log((1+df$FemaleEnriched)/(1+df$MaleEnriched))

# Store MGI into tmp phylobject
tmp@sam_data$MGI <- df$MGI


# Computing Alpha diversity using Observed diversity
alpha <- plot_richness(df_alpha, x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
df <- alpha[["data"]]
df$MGI <- tmp@sam_data$MGI


Plot_MGI_Perio <- ggscatter(df, x = "MGI", y = "value",
                    add = "reg.line",                                       # Add regression line
                    conf.int = TRUE,                                        # Add confidence interval
                    fullrange = TRUE,
                    palette = c("#F88B9D", "#8ECEFD"),
                    alpha = 0.6,
                    shape = 21,
                    size = 3,
                    fill = "Sex",
                    color = "Sex",
                    add.params = list(color = "#252F63", fill = "gray90", size = .5))+
                    labs(x="MG-DI (Log)", y="Observed ASV diversity")+
                    MyTheme+
                    stat_cor(method = "pearson", label.x = 0.1, label.y = 3.1)
# Store my enriched genera plot 
MLresults[["Plot_MGI_Perio"]] <- Plot_MGI_Perio

# Verify AUCROC of MG-DI with ML
library(caret)
library(MLeval)

foo <- na.omit(df[, c("Sex", "MGI")])
foo$Sex <- factor(foo$Sex)

# Setting seed for reproducibility
set.seed(13122019)

# Setting the ML
cntrl <- trainControl(method = "repeatedcv",
                      number = 10, repeats=3,    # 3 repeats of 10-fold cross-validation, it will perform 10-fold cross-validation on the training data 3 times
                      summaryFunction=prSummary,
                      classProbs=T,
                      savePredictions = T,
                      verboseIter = TRUE)

# Running Random Forest
rf <- train(Sex ~., 
            data = foo, 
            method = "rf",
            metric = "AUC",
            trControl = cntrl)

## Analyze the ML models with  MLeval
res <- evalm(list(rf), showplots = FALSE, title = "Periodontitis",
             gnames = c("RF"),
             cols = c("#999999"))

# Store my ML resuls
MLresults[["PlotPerio"]] <- res$roc$data
MLresults[["MLPerio"]] <- data.frame(res$stdres$RF, "Var" = rownames(res$stdres$RF))







foo <- na.omit(df[, c("Sex", "MGI", "AgeYrs")])
foo$Sex <- factor(foo$Sex)

# Setting seed for reproducibility
set.seed(13122019)

# Setting the ML
cntrl <- trainControl(method = "repeatedcv",
                      number = 10, repeats=3,    # 3 repeats of 10-fold cross-validation, it will perform 10-fold cross-validation on the training data 3 times
                      summaryFunction=prSummary,
                      classProbs=T,
                      savePredictions = T,
                      verboseIter = TRUE)

# Running Random Forest
rf <- train(Sex ~., 
            data = foo, 
            method = "rf",
            metric = "AUC",
            trControl = cntrl)

## Analyze the ML models with  MLeval
res <- evalm(list(rf), showplots = FALSE, title = "Periodontitis and age",
             gnames = c("RF"),
             cols = c("#999999"))


# Store my ML resuls
# Store my ML resuls
MLresults[["PlotPerioAge"]] <- res$roc$data
MLresults[["MLPerioAge"]] <- data.frame(res$stdres$RF, "Var" = rownames(res$stdres$RF))



