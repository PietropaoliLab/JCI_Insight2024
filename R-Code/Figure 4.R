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


load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
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
                            labs(x="", y="Effect size estimated mean", title = "Periodontitis")+
                            coord_flip(ylim = c(-0.010, 0.010))+
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
  labs(x="MGI (Log)", y="Observed ASV diversity", title = "Periodontitis")+
  MyTheme+
  scale_y_continuous(breaks = seq(0, 60, by=10))+
  coord_cartesian(ylim = c(0, 60))+
  stat_cor(method = "pearson", label.x = 0.1, label.y = 3.1)

# Store my enriched genera plot 
MLresults[["Plot_MGI_Perio"]] <- Plot_MGI_Perio


MGI_Age_Perio <- ggscatter(df, x = "AgeYrs",y = "MGI",
                           add = "reg.line",                                       # Add regression line
                           conf.int = TRUE,                                        # Add confidence interval
                           fullrange = TRUE,
                           palette = c("#F88B9D", "#8ECEFD"),
                           alpha = 0.6,
                           shape = 21,
                           size = 3,
                           fill = "Sex",
                           color = "Sex")+
                          labs(x="Age", y="MGI (Log)", title = "Periodontitis")+
                          MyTheme+
                          #scale_y_continuous(breaks = seq(0, 60, by=10))+
                          #coord_cartesian(ylim = c(0, 60))+
                          stat_cor(aes(color = Sex), method = "pearson")

MLresults[["MGI_Age_Perio"]] <- MGI_Age_Perio




MGI_Perio <- ggplot(df, aes(x=Sex, y=MGI, color = Sex))+             
                    geom_violin(aes(fill=Sex), alpha=.3)+
                    geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=3, stroke = .7)+
                    geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
                    labs(title = "Periodontitis", y = "MGI", x="")+
                    #scale_y_continuous(breaks = c(20, 40, 60, 80))+
                    #coord_fixed(ylim = c(10,90))+
                    scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
                    scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
                    MyTheme+
                    theme(aspect.ratio = 2/1)+
                    stat_compare_means(comparisons = list(c("female", "male")), method = "wilcox.test", label = "p.signif", family="Avenir")

MLresults[["MGI_Perio"]] <- MGI_Perio


# Verify AUCROC of MGI with ML
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
MLresults[["PlotPerioAge"]] <- res$roc$data
MLresults[["MLPerioAge"]] <- data.frame(res$stdres$RF, "Var" = rownames(res$stdres$RF))




################################################################################

#                  H E A L T H Y 

################################################################################
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# my results df
results <- data.frame()

# Define my subseto of interest
tmp <- subset_samples(psm_healthy, Site == "Subgingival" & Smoking == "No")

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
PlotHealthyEnrichedGenera <- ggplot(results, aes(x=reorder(feature, ef_diff_mean), y = ef_diff_mean, fill = enrich_group))+
                                    geom_hline(yintercept = 0, size = .1)+
                                    geom_col(width = .75)+
                                    scale_fill_manual(values = MyPalette)+
                                    labs(x="", y="Effect size estimated mean", title = "Healthy")+
                                    coord_flip(ylim = c(-0.010, 0.010))+
                                    MyTheme
PlotHealthyEnrichedGenera

# Store my enriched genera plot 
MLresults[["PlotHealthyEnrichedGenera"]] <- PlotHealthyEnrichedGenera

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

# Computing Alpha diversity using Shannon
alpha <- plot_richness(df_alpha, x="Sex", measures="Observed", color = "Sex")

# Keep my data into a new dataset
df <- alpha[["data"]]
df$MGI <- tmp@sam_data$MGI



Plot_MGI_Healthy <- ggscatter(df, x = "MGI", y = "value",
                          add = "reg.line",                                       # Add regression line
                          conf.int = TRUE,                                        # Add confidence interval
                          fullrange = TRUE,
                          shape = 21,
                          palette = c("#F88B9D", "#8ECEFD"),
                          alpha = 0.6,
                          size = 3,
                          fill = "Sex",
                          color = "Sex",
                          add.params = list(color = "#252F63", fill = "gray90", size = .5))+
  labs(x="MGI (Log)", y="Observed ASV diversity", title = "Healthy")+
  MyTheme+
  coord_cartesian(ylim = c(0, 60))+
  scale_y_continuous(breaks = seq(0, 60, by=10))+
  stat_cor(method = "pearson", label.x = 0.001, label.y = 2.6)

# Store my enriched genera plot 
MLresults[["Plot_MGI_Healthy"]] <- Plot_MGI_Healthy


MGI_Age_Healthy <- ggscatter(df, x = "AgeYrs", y = "MGI",
                           add = "reg.line",                                       # Add regression line
                           conf.int = TRUE,                                        # Add confidence interval
                           fullrange = TRUE,
                           palette = c("#F88B9D", "#8ECEFD"),
                           alpha = 0.6,
                           shape = 21,
                           size = 3,
                           fill = "Sex",
                           color = "Sex")+
                           labs(x="Age", y="MGI (Log)", title = "Healthy")+
                            MyTheme+
                            #scale_y_continuous(breaks = seq(0, 60, by=10))+
                            #coord_cartesian(ylim = c(0, 60))+
                            stat_cor(aes(color = Sex), method = "pearson")

MLresults[["MGI_Age_Healthy"]] <- MGI_Age_Healthy



MGI_Healthy <- ggplot(df, aes(x=Sex, y=MGI, color = Sex))+             
                      geom_violin(aes(fill=Sex), alpha=.3)+
                      geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=3, stroke = .7)+
                      geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
                      labs(title = "Healthy", y = "MGI", x="")+
                      #scale_y_continuous(breaks = c(20, 40, 60, 80))+
                      #coord_fixed(ylim = c(10,90))+
                      scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
                      scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
                      MyTheme+
                      theme(aspect.ratio = 2/1)+
                      stat_compare_means(comparisons = list(c("female", "male")), method = "wilcox.test", label = "p.signif", family="Avenir")

MLresults[["MGI_Healthy"]] <- MGI_Healthy


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
res <- evalm(list(rf), showplots = FALSE, title = "Healthy",
             gnames = c("RF"),
             cols = c("#999999"))


# Store my ML resuls
MLresults[["PlotHealthy"]] <- res$roc$data
MLresults[["MLHealthy"]] <- data.frame(res$stdres$RF, "Var" = rownames(res$stdres$RF))




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
res <- evalm(list(rf), showplots = FALSE, title = "Healthy and age",
             gnames = c("RF"),
             cols = c("#999999"))



# Store my ML resuls
MLresults[["PlotHealthyAge"]] <- res$roc$data
MLresults[["MLHealthyAge"]] <- data.frame(res$stdres$RF, "Var" = rownames(res$stdres$RF))

# Cleaning
# rm(alpha, cntrl, df, df_alpha, FemaleGenus, foo, Genera, 
#   MaleGenus, MGI, MyPalette, 
#   psm_healthy, psm_perio, res, results, rf, tmp)


# Plotting ROCs
# periodontitis
MGIalone <- paste0("AUC-ROC:", subset(MLresults$MLPerio, Var == "AUC-ROC")$Score)
MGIwAge <- paste0("AUC-ROC:", subset(MLresults$MLPerioAge, Var == "AUC-ROC")$Score)

ROC_Perio <- ggplot()+
                  geom_abline(intercept = 0, slope = 1, size = 0.25, color = "gray 60") +
            
                  geom_line(data = MLresults$PlotPerio, aes(x = FPR, y= SENS, color = "MGI alone"),  size = 2)+ # Only Periodontitis
                  geom_line(data = MLresults$PlotPerioAge, aes(x = FPR, y= SENS, color = "MGI+age"),  size = 2)+ # Only Periodontitis and Age
                  # Set manual colors and legend
                  scale_color_manual(name='ML', values=c('MGI alone'='#ff8288', 'MGI+age'='#b3494e'))+
                  
                  # MGI alone
                  annotate("text", x = 1, y = 0, vjust=0, hjust = 1, label = MGIalone, color = '#ff8288', family = "Avenir", size = 3)+
                  # MGI with age
                  annotate("text", x = 1, y = 0.10, vjust=0, hjust = 1, label = MGIwAge, color = '#b3494e', family = "Avenir", size = 3)+
                  
                  labs(x="False positive rate", y="True positive rate", title = "Periodontitis")+
                  MyTheme+
                  theme(legend.position = "bottom")


# Healthy
MGIalone <- paste0("AUC-ROC:", subset(MLresults$MLHealthy, Var == "AUC-ROC")$Score)
MGIwAge <- paste0("AUC-ROC:", subset(MLresults$MLHealthyAge, Var == "AUC-ROC")$Score)

ROC_Healthy <- ggplot()+
                    geom_abline(intercept = 0, slope = 1, size = 0.25, color = "gray 60") +
                    
                    geom_line(data = MLresults$PlotHealthy, aes(x = FPR, y= SENS, color = "MGI alone"),  size = 2)+ # Only Periodontitis
                    geom_line(data = MLresults$PlotHealthyAge, aes(x = FPR, y= SENS, color = "MGI+age"),  size = 2)+ # Only Periodontitis and Age
                    # Set manual colors and legend
                    scale_color_manual(name='ML', values=c('MGI alone'='#04bf55', 'MGI+age'='#017332'))+
                    
                    # MGI alone
                    annotate("text", x = 1, y = 0, vjust=0, hjust = 1, label = MGIalone, color = '#04bf55', family = "Avenir", size = 3)+
                    # MGI with age
                    annotate("text", x = 1, y = 0.10, vjust=0, hjust = 1, label = MGIwAge, color = '#017332', family = "Avenir", size = 3)+
                    
                    labs(x="False positive rate", y="True positive rate", title = "Healthy")+
                    MyTheme+
                    theme(legend.position = "bottom")






library(patchwork)
# Combine my genera
Genera <- (MLresults$PlotHealthyEnrichedGenera | MGI_Healthy) / (MLresults$PlotPerioEnrichedGenera | MGI_Perio)
Genera

# Combine my MGI
MGI <- (Plot_MGI_Healthy | Plot_MGI_Perio)
MGI

# ROC
ROC <- (ROC_Healthy | ROC_Perio)
ROC


# Regression MGI Age
RegMGI_AGE <- (MGI_Age_Healthy | MGI_Age_Perio)
RegMGI_AGE




ggsave(plot = Genera, 
       device = cairo_pdf,
       width = 8, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Genera and Violin Figure 4.pdf")


ggsave(plot = MGI, 
       device = cairo_pdf,
       width = 4, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/MGI Figure 4.pdf")



ggsave(plot = ROC, 
       device = cairo_pdf,
       width = 4, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/ROC Figure 4.pdf")



ggsave(plot = RegMGI_AGE, 
       device = cairo_pdf,
       width = 4, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/RegMGI_AGE Figure 4.pdf")





EllipseH <- ggscatter(Plot_MGI_Healthy$data, 
          x = "MGI", y = "value",
          #add = "reg.line",                                       # Add regression line
          conf.int = FALSE,                                        # Add confidence interval
          fullrange = TRUE,
          ellipse = TRUE,
          shape = 21,
          palette = c("#F88B9D", "#8ECEFD"),
          alpha = 0.6,
          size = 3,
          fill = "Sex",
          color = "Sex")+
          labs(x="MGI (Log)", y="Observed ASV diversity", title = "Healthy")+
          MyTheme+
          coord_cartesian(ylim = c(0, 65))+
          scale_y_continuous(breaks = seq(0, 60, by=10))+
          stat_cor(aes(color = Sex), method = "pearson", label.x = 0, label.y = c(5, 10))


EllipseP <- ggscatter(Plot_MGI_Perio$data, 
          x = "MGI", y = "value",
          #add = "reg.line",                                       # Add regression line
          conf.int = FALSE,                                        # Add confidence interval
          fullrange = TRUE,
          ellipse = TRUE,
          shape = 21,
          palette = c("#F88B9D", "#8ECEFD"),
          alpha = 0.6,
          size = 3,
          fill = "Sex",
          color = "Sex")+
          labs(x="MGI (Log)", y="Observed ASV diversity", title = "Periodontitis")+
          MyTheme+
          coord_cartesian(ylim = c(0, 65))+
          scale_y_continuous(breaks = seq(0, 60, by=10))+
          stat_cor(aes(color = Sex), method = "pearson", label.x = 0.1, label.y = c(5, 10))


Ellipse <- EllipseH | EllipseP

ggsave(plot = Ellipse, 
       device = cairo_pdf,
       width = 4, 
       #height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Ellipse Figure 4.pdf")







