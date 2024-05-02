library(phyloseq) 
library(microbiome)
library(caret)
library(MLeval)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(microbiomeMarker)
library(extrafont)
library(showtext)
loadfonts(device = "pdf")


# OSX
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")




######################################################
#
#         WITHIN STUDIES SAMPLE SIZE ANALYSIS
#
######################################################

healthy <- psm_healthy@sam_data
perio <- psm_perio@sam_data

heathy_perio <- rbind(healthy, perio)

Freq <- ddply(heathy_perio, c('BioProject','Site', 'Periodontitis'), summarise,
              Freqs = sum(!is.na(BioSample)))

Freq$Label <- paste0(Freq$BioProject, " - ", Freq$Periodontitis)
Freq$LongLabel <- paste0(Freq$BioProject, " - ", Freq$Periodontitis, " - ", Freq$Site)

# sorting
attach(Freq)
Freq <- Freq[order(Site, Periodontitis, Freqs),]
detach(Freq)

Freq$Ordering <- seq(1:nrow(Freq))






######################################################
#
#         WITHIN STUDIES RANDOM FOREST
#
######################################################

# Staring cycle for Periodontitis
bioproject <- unique(psm_perio@sam_data$BioProject)
sc.list <- list()
aucroc <- data.frame()

for (b in bioproject) {

tmp <- psm_perio
tmp <- tax_glom(tmp, taxrank = "Genus")
tmp <- subset_samples(tmp, BioProject == b)
site <- unique(tmp@sam_data$Site)


# Clycle for periodontitis
for (s in site) {
  foo <- phyloseq::subset_samples(tmp, Site == s)
  
  
  # Melting Phyloseq Objects
  foo <- psmelt(foo)
  
  # My vars for training ML
  myvars<-c("Run","Sex", "Genus", "Abundance") 
  
  
  # Define the DF of interest providing NA cleaning
  df<-foo[myvars]
  df<-na.omit(df)
  
  # Generate abundance table for ML
  df$AbundanceLog10 <- log10(df$Abundance + 1)
  abn <- df[c("Run", "Genus", "AbundanceLog10")]
  
  # Fixing coloum names
  abn$Genus <- gsub(" ", "_", abn$Genus)
  abn$Genus <- gsub("-", "_", abn$Genus)
  abn$Genus <- gsub("\\[", "", abn$Genus)
  abn$Genus <- gsub("\\]", "", abn$Genus)
  abn$Genus <- gsub("\\.", "", abn$Genus)
  
  # Generate abundance wide for further ML steps
  abn_wide <- dcast(abn, Run ~ Genus, value.var="AbundanceLog10")
  
  # Subsetting variable of interest
  df <- df[c("Run", "Sex")]
  df <- unique(df)
  
  # Merging Metadata and abundance wide table
  df <- merge(df, abn_wide, by = "Run")
  rownames(df) <- df$Run
  
  # Setting factors
  df$Sex <- factor(df$Sex)
  #df$Smoking <- factor(df$Smoking)
  #df$Site <- factor(df$Site)
  
  # Purging the dataset
  df <- na.omit(df)
  
  # Remove the "Run" coloum -- IDs
  df <- df[-1]
  
  # Setting seed for reproducibility
  set.seed(101)
  
  # Setting the ML
  cntrl <- trainControl(method = "repeatedcv",
                        number = 10, repeats=3,    # 3 repeats of 10-fold cross-validation, it will perform 10-fold cross-validation on the training data 3 times
                        summaryFunction=prSummary,
                        classProbs=T,
                        savePredictions = T,
                        verboseIter = TRUE)
  
  
 
  # Running Random Forest
  rf <- train(Sex ~., 
              data = df, 
              method = "rf",
              metric = "AUC",
              trControl = cntrl)
  
  
  
  ## Analyze the ML models with  MLeval
  res <- evalm(list(rf), showplots = FALSE, title = paste0(b, " perio-", s),
               gnames = c("RF"),
               cols = c("#999999"))
  
  x <- as.data.frame(res$optres)[13,]
  x$BioProject <- b
  x$Site <- s
  x$Periodontitis <- 'Periodontitis'
  aucroc <- rbind(aucroc, x)
  
  sc.list[[paste0(b,"perio-", s)]] <- res
  
}

}



# Staring cycle for Healthy
bioproject <- unique(psm_healthy@sam_data$BioProject)


for (b in bioproject) {
  
  tmp <- psm_healthy
  tmp <- tax_glom(tmp, taxrank = "Genus")
  tmp <- subset_samples(tmp, BioProject == b)
  site <- unique(tmp@sam_data$Site)
  
  
  # Clycle for healthy
  for (s in site) {
    foo <- phyloseq::subset_samples(tmp, Site == s)
    
    
    # Melting Phyloseq Objects
    foo <- psmelt(foo)
    
    # My vars for training ML
    myvars<-c("Run","Sex", "Genus", "Abundance") 
    
    
    # Define the DF of interest providing NA cleaning
    df<-foo[myvars]
    df<-na.omit(df)
    
    # Generate abundance table for ML
    df$AbundanceLog10 <- log10(df$Abundance + 1)
    abn <- df[c("Run", "Genus", "AbundanceLog10")]
    
    # Fixing coloum names
    abn$Genus <- gsub(" ", "_", abn$Genus)
    abn$Genus <- gsub("-", "_", abn$Genus)
    abn$Genus <- gsub("\\[", "", abn$Genus)
    abn$Genus <- gsub("\\]", "", abn$Genus)
    abn$Genus <- gsub("\\.", "", abn$Genus)
    
    # Generate abundance wide for further ML steps
    abn_wide <- dcast(abn, Run ~ Genus, value.var="AbundanceLog10")
    
    # Subsetting variable of interest
    df <- df[c("Run", "Sex")]
    df <- unique(df)
    
    # Merging Metadata and abundance wide table
    df <- merge(df, abn_wide, by = "Run")
    rownames(df) <- df$Run
    
    # Setting factors
    df$Sex <- factor(df$Sex)
    #df$Smoking <- factor(df$Smoking)
    #df$Site <- factor(df$Site)
    
    # Purging the dataset
    df <- na.omit(df)
    
    # Remove the "Run" coloum -- IDs
    df <- df[-1]

    # Setting seed for reproducibility
    set.seed(101)
    
    # Setting the ML
    cntrl <- trainControl(method = "repeatedcv",
                          number = 9, repeats=3,    # 3 repeats of 10-fold cross-validation, it will perform 10-fold cross-validation on the training data 3 times
                          summaryFunction=prSummary,
                          classProbs=T,
                          savePredictions = T,
                          verboseIter = TRUE)
    
    
    
    # Running Random Forest
    rf <- train(Sex ~., 
                data = df, 
                method = "rf",
                metric = "AUC",
                trControl = cntrl)
    
    
    
    ## Analyze the ML models with  MLeval
    res <- evalm(list(rf), showplots = FALSE, title = paste0(b, " healthy-", s),
                 gnames = c("RF"),
                 cols = c("#999999"))
    
    x <- as.data.frame(res$optres)[13,]
    x$BioProject <- b
    x$Site <- s
    x$Periodontitis <- 'Healthy'
    aucroc <- rbind(aucroc, x)
    
    sc.list[[paste0(b,"healthy-", s)]] <- res
    
  }
  
}






######################################################
#
#    EVALUATING GENERA DIFFERENCES WITHIN STUDIES
#           USING welch.test AND FDR
#
######################################################

# Setting my variables
Sites <- c("Saliva", "Plaque", "Subgingival")
res <- data.frame()

# # ================================
# # Clycle for periodontitis
# # ================================
# Within Sites
for (s in Sites) {
  foo <- subset_samples(psm_perio, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  
  # Within dataset
  for (d in datasets){
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Print some info on screen
    cat(paste0("Analyzing: ", d, "-", s, "\n"))
    
    # Operate my analysis
    tmp <- run_simple_stat(ps,
                           taxa_rank = "Genus",
                           transform = "log10p",
                           norm = "TSS",
                           method = "welch.test",
                           # p_adjust = "fdr",
                           pvalue_cutoff = 1, # Incorporate all Phyla
                           group = "Sex")
    # My results dataset
    tmp_res <- tmp@marker_table
    
    # Adding FDR to p-value
    tmp_res$padj <- p.adjust(tmp_res$pvalue, method = "fdr")
    
    # Adding type of test used
    tmp_res$Stat_method <- tmp@diff_method
    
    # Adding groups
    tmp_res$BioProject <- d
    tmp_res$Site <- s
    tmp_res$Periodontitis <- "Periodontitis"
    
    # Storing results
    res <- rbind(res, tmp_res)
    rm(tmp_res, tmp)
    
  }
}
rm(foo, ps, d, s)



# # ================================
# # Clycle for Healthy
# # ================================
# Within Sites
for (s in Sites) {
  foo <- subset_samples(psm_healthy, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  
  # Within dataset
  for (d in datasets){
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Print some info on screen
    cat(paste0("Analyzing: ", d, "-", s, "\n"))
    
    # Operate my analysis
    tmp <- run_simple_stat(ps,
                           taxa_rank = "Genus",
                           transform = "log10p",
                           norm = "TSS",
                           method = "welch.test",
                           # p_adjust = "fdr",
                           pvalue_cutoff = 1, # Incorporate all Phyla
                           group = "Sex")
    # My results dataset
    tmp_res <- tmp@marker_table
    
    # Adding FDR to p-value
    tmp_res$padj <- p.adjust(tmp_res$pvalue, method = "fdr")
    
    # Adding type of test used
    tmp_res$Stat_method <- tmp@diff_method
    
    # Adding groups
    tmp_res$BioProject <- d
    tmp_res$Site <- s
    tmp_res$Periodontitis <- "Healthy"
    
    # Storing results
    res <- rbind(res, tmp_res)
    rm(tmp_res, tmp)
    
  }
}
rm(foo, ps, d, s, Sites, datasets) # some cleaning





######################################################
#
#      WITHIN STUDIES MICROGENDEROME INDEX
#
######################################################

MGI <- data.frame()

bioproject <- unique(psm_perio@sam_data$BioProject)

for (b in bioproject) {
#  b="PRJNA774981"
  tmp <- psm_perio
  tmp <- tax_glom(tmp, taxrank = "Genus")
  tmp <- subset_samples(tmp, BioProject == b)
  site <- unique(tmp@sam_data$Site)
  
  cat(b, ":")
  # Clycle for Perio
  for (s in site) {
   # s="Saliva"
    # Take significant Genus from KW results
    MyGenus <- subset(KW, BioProject == b & Site == s & Periodontitis == "Periodontitis")
    GenusFemale <- subset(MyGenus, enrich_group == "female")$Genus
    GenusMale <- subset(MyGenus, enrich_group == "male")$Genus
    
    # Subsetting dataset according to BioProject and Sampling Site
    print(table(tmp@sam_data$BioProject, tmp@sam_data$Site))
    
    # Relative abundances
    tmp <- microbiome::transform(tmp, "compositional")
    
    # Computing Enriched genera in Female
    if (length(GenusFemale) <1 ) { # If no Genus were present TaxEnrich == 1
      TaxEnrichFemale <- 1
      
    } else {
      TaxEnrichFemale <- subset_taxa(tmp, Genus %in% GenusFemale)
      TaxEnrichFemale <- psmelt(TaxEnrichFemale)
      TaxEnrichFemale <- (sum(TaxEnrichFemale$Abundance, na.rm = TRUE)+1)
    }
    
    
    # Computing Enriched genera in Male
    if (length(GenusMale) <1 ) { # If no Genus were present TaxEnrich == 1
      TaxEnrichMale <- 1
      
    } else {
      TaxEnrichMale <- subset_taxa(tmp, Genus %in% GenusMale)
      TaxEnrichMale <- psmelt(TaxEnrichMale)
      TaxEnrichMale <- (sum(TaxEnrichMale$Abundance, na.rm = TRUE)+1)
    }
    
    # Generating dataset of results
    foo <- data.frame(BioProject = b,
                      Site = s,
                      Periodontitis = "Periodontitis",
                      TaxEnrichFemale = TaxEnrichFemale,
                      TaxEnrichMale = TaxEnrichMale)
    
    foo$MGI_raw = foo$TaxEnrichFemale/foo$TaxEnrichMale
    foo$MGI_Log = log(foo$MGI_raw)
    foo$MGI_Log10 = log10(foo$MGI_raw)
    
    
    print(foo)
    # Store the results
    MGI <- rbind(MGI, foo)
  }
}





# MGI for healthy
bioproject <- unique(psm_healthy@sam_data$BioProject)

for (b in bioproject) {
  # b="PRJEB6047"
  tmp <- psm_healthy
  tmp <- tax_glom(tmp, taxrank = "Genus")
  tmp <- subset_samples(tmp, BioProject == b)
  site <- unique(tmp@sam_data$Site)
  
  cat(b, ":")
  # Clycle for Healthy
  for (s in site) {
    # s="Plaque"
    # Take significant Genus from KW results
    MyGenus <- subset(KW, BioProject == b & Site == s & Periodontitis == "Healthy")
    GenusFemale <- subset(MyGenus, enrich_group == "female")$Genus
    GenusMale <- subset(MyGenus, enrich_group == "male")$Genus
    
    # Subsetting dataset according to BioProject and Sampling Site
    print(table(tmp@sam_data$BioProject, tmp@sam_data$Site))
    
    # Relative abundances
    tmp <- microbiome::transform(tmp, "compositional")
    
    # Computing Enriched genera in Female
    if (length(GenusFemale) <1 ) { # If no Genus were present TaxEnrich == 1
      TaxEnrichFemale <- 1
      
    } else {
      TaxEnrichFemale <- subset_taxa(tmp, Genus %in% GenusFemale)
      TaxEnrichFemale <- psmelt(TaxEnrichFemale)
      TaxEnrichFemale <- (sum(TaxEnrichFemale$Abundance, na.rm = TRUE)+1)
    }
    
    
    # Computing Enriched genera in Male
    if (length(GenusMale) <1 ) { # If no Genus were present TaxEnrich == 1
      TaxEnrichMale <- 1
      
    } else {
      TaxEnrichMale <- subset_taxa(tmp, Genus %in% GenusMale)
      TaxEnrichMale <- psmelt(TaxEnrichMale)
      TaxEnrichMale <- (sum(TaxEnrichMale$Abundance, na.rm = TRUE)+1)
    }
    
    # Generating dataset of results
    foo <- data.frame(BioProject = b,
                      Site = s,
                      Periodontitis = "Healthy",
                      TaxEnrichFemale = TaxEnrichFemale,
                      TaxEnrichMale = TaxEnrichMale)
    
    foo$MGI_raw = foo$TaxEnrichFemale/foo$TaxEnrichMale
    foo$MGI_Log = log(foo$MGI_raw)
    foo$MGI_Log10 = log10(foo$MGI_raw)
    
    
    print(foo)
    # Store the results
    MGI <- rbind(MGI, foo)
  }
}







######################################################
#
#             COMPOSING PLOTS
#
######################################################

library(tidytext)
library(ggh4x)

# -------------------#
## Sample size plot
# -------------------#
# Ordering my sites
Freq$Site <- factor(Freq$Site, levels = c("Saliva", "Plaque", "Subgingival"))

p1 <- ggplot(Freq, aes(x=reorder(LongLabel, Ordering), y=Freqs, fill=Periodontitis)) + 
  geom_hline(yintercept = c(20, 60, 100), size = .2, linetype =2, color = 'gray65')+
  geom_bar(position=position_dodge(.2), stat="identity", width = .65) + 
  facet_grid(Site~.,  scales = "free", switch = "y") + 
  labs(x='', y='Samples')+
  theme_classic() +
  scale_fill_manual(values = c("#41ab74", "#eb5c60")) +
  coord_flip(ylim = c(0, 102))+
  scale_y_continuous(breaks = c(20, 60,  100), expand=c(0,0))+
  scale_x_reordered()+
  theme_classic(base_size = 11)+
  theme(legend.title=element_blank(), legend.position="none", 
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
        panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),   
        strip.text.y = element_text(angle=0, size = 13),
        strip.background.y=element_part_rect(side = "b", colour = "black", fill = NA),
        strip.placement = "outside",
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        aspect.ratio = .8/1,
        text = element_text(family = "Avenir", color = "black"),
        axis.text = element_text(colour = "black", size = 10))


# -------------------#
## AUCROC plot
# -------------------#
# Operate on AUCROC dataset
aucroc$LongLabel <- paste0(aucroc$BioProject, " - ", aucroc$Periodontitis, " - ", aucroc$Site)
aucroc$Label <- paste0(aucroc$BioProject, " - ", aucroc$Periodontitis)
aucroc <- merge(aucroc, Freq[,c('LongLabel', 'Ordering')])
aucroc$Site <- factor(aucroc$Site, levels = c("Saliva", "Plaque", "Subgingival"))

p2 <- ggplot(aucroc, aes(x=reorder(LongLabel, Ordering), y=RF.Score)) + 
  geom_hline(yintercept = c(.25, .50, .75), size = .2, linetype =2, color = 'gray65')+

  geom_segment(aes(x=reorder(LongLabel, Ordering), 
                   xend=reorder(LongLabel, Ordering), 
                   y=0, 
                   yend=RF.Score), size = 1.2, color = 'gray85') +
  
  geom_point(aes(x=reorder(LongLabel, Ordering), 
                 y=RF.Score), size=4.2, shape = 21, color = 'gray85', fill ='gray85') +
  geom_point(aes(x=reorder(LongLabel, Ordering), 
                 y=RF.Score, fill=Periodontitis), size = 2.8, shape = 21, color = 'white')+
    facet_grid(Site~.,  scales = "free") + 
    labs(x='', y='AUC-ROC')+
  scale_y_continuous(breaks = c(0, .50,  1), expand=c(0,0))+
   # scale_x_reordered()+
  scale_fill_manual(values = c("#41ab74", "#eb5c60")) +
  coord_flip(ylim = c(0, 1.02))  +
  theme_classic(base_size = 11)+
  theme(legend.title=element_blank(), legend.position="none", 
        # axis.ticks.x=element_blank(), 
        # axis.ticks.y=element_blank(), 
        axis.text.y = element_blank(), # hide axis text
        panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),   
        strip.text.y = element_blank(), # hide labels
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        aspect.ratio = .8/1,
        text = element_text(family = "Avenir", color = "black"),
        axis.text = element_text(colour = "black", size = 10))



# ----------------------#
## welch.test plot
# ----------------------#

# Operate on welch.test dataset
res$LongLabel <- paste0(res$BioProject, " - ", res$Periodontitis, " - ", res$Site)
res$Label <- paste0(res$BioProject, " - ", res$Periodontitis)
res <- merge(res, Freq[,c('LongLabel', 'Ordering')], all.y = TRUE)

# Compute Significant Genera
welch <- ddply(res, c("LongLabel", "Label", "Ordering", "enrich_group", "Site"), summarise,
               N = sum(padj <0.05))

welch$Sex <- ifelse(welch$enrich_group == "female", "F", "M")

# Set female on left side
welch$N <- ifelse(welch$Sex == "F", - abs(welch$N), welch$N)

welch$Site <- factor(welch$Site, levels = c("Saliva", "Plaque", "Subgingival"))



p3 <- ggplot(welch, aes(group = Sex)) + 
  geom_hline(yintercept = c(0), size = .4, linetype =1, color = 'black')+
  geom_hline(yintercept = c(-10, -5, 5, 10), size = .2, linetype =2, color = 'gray65')+
  
  # Line of Lolliplot 
  geom_linerange(aes(y = N, 
                     ymin = 0,
                     ymax = N,
                     x = reorder(LongLabel, Ordering), color = Sex),
                     position = position_dodge(.0), size = 1.2)+ 
  # External dot
  geom_point(aes(x = reorder(LongLabel, Ordering), y = N, color = Sex), 
                position = position_dodge(.0),
                size= 2.8, shape = 16) +
  
    scale_color_manual(values = c('#f2acda', '#6c9bd9'))+

  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), labels = c("10", "5", "0", "5", "10"))+
  
  facet_grid(Site~.,  scales = "free", drop=TRUE) + 
  labs(x='', y='Genus enriched')+
  
  theme_classic(base_size = 11)+
  coord_flip(ylim = c(-11, 13.5))  +
  theme(legend.title=element_blank(), legend.position="none", 
        axis.text.y = element_blank(), # hide axis text
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
        panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),   
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        aspect.ratio = .8/1,
        text = element_text(family = "Avenir", color = "black"),
        axis.text = element_text(colour = "black", size = 10))  




# ----------------------#
## Microgenderome Index plot
# ----------------------#
# Operate on MGI dataset
MGI$LongLabel <- paste0(MGI$BioProject, " - ", MGI$Periodontitis, " - ", MGI$Site)
MGI$Label <- paste0(MGI$BioProject, " - ", MGI$Periodontitis)
MGI <- merge(MGI, Freq[,c('LongLabel', 'Ordering')], all.y = TRUE)

p4 <- ggplot(MGI, fill = site) + 
  geom_hline(yintercept = c(-1.5, 1.5), size = .2, linetype =2, color = 'gray65')+
  geom_hline(yintercept = 0, size = .2, linetype =1, color = 'gray65')+
  
    geom_segment(aes(x=reorder_within(Label, Ordering, Site), 
                   xend=reorder_within(Label, Ordering, Site), 
                   y=0, 
                   yend=MGI_Log10), size = 1.2, color = 'gray85') +
  
  geom_point(aes(x=reorder_within(Label, Ordering, Site), 
                 y=MGI_Log10), size=4.2, shape = 21, color = 'gray85', fill ='gray85') +
  geom_point(aes(x=reorder_within(Label, Ordering, Site), 
                 y=MGI_Log10, fill=Site), size = 2.8, shape = 21, color = 'white')+
  
  
  scale_y_continuous(breaks = c(-1, -0.5, 0, .5, 1), expand=c(0,0))+
  scale_fill_brewer(palette="Set2") +
  
  facet_grid(Site~.,  scales = "free", drop=TRUE) + 
  labs(x='', y='MGI index (Log10)')+
  theme_classic() +
  scale_x_reordered()+
  coord_flip(ylim = c(-1.55, 1.55))  +
  theme(legend.title=element_blank(), legend.position="none", 
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
        panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),   
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        aspect.ratio = 1/1)  







library(patchwork)
MetaAnalysis <- p1 + p2 + p3

MetaAnalysis
ggsave(plot = MetaAnalysis, 
       device = cairo_pdf,
       # width = 5.5, 
       # height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/MetaAnalysis.pdf")


MetaAnalysis_SampleSize <- p1
MetaAnalysis_AUC <- p2
MetaAnalysis_GenusEnriched <- p3






