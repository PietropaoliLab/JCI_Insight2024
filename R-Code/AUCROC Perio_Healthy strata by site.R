library(phyloseq) 
library(microbiome)
library(caret)
library(MLeval)
library(reshape2)
library(ggplot2)
library(extrafont)
library(showtext)
# library(ggprism)


loadfonts(device = "pdf")


# Load my data
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# # windows
# load("G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
# load("G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank

  site <- c("Saliva", "Plaque", "Subgingival")
  sc.list <- list()

  tmp <- psm_perio
  tmp <- tax_glom(tmp, taxrank = "Genus")

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


# Training ML Genus and Smoking habits
# Running GBM here
gbm <- train(Sex ~., 
             data = df, 
             method = "gbm",
             metric = "AUC",
             bag.fraction = 0.7, # bag.fraction (Subsampling fraction) - the fraction of the training set observations randomly selected to propose the next tree in the expansion. In this case, it adopts stochastic gradient boosting strategy. By default, it is 0.5. That is half of the training sample at each iteration. You can use fraction greater than 0.5 if training sample is small.
             trControl = cntrl)

# Running Random Forest
rf <- train(Sex ~., 
            data = df, 
            method = "rf",
            metric = "AUC",
            trControl = cntrl)


# ========================================= #
#       Clycle for analyzing both 
#         Periodontitis & Healthy
#
# ========================================= #
## Analyze the ML models with  MLeval
res <- evalm(list(gbm, rf), showplots = TRUE, title = paste0("Periodontitis - ", s),
             gnames = c("GBM", "RF"),
             cols = c("#ca0020", "#999999"))


sc.list[[paste0("Periodontitis - ", s)]] <- res

}

 
  
  tmp <- psm_healthy
  tmp <- tax_glom(tmp, taxrank = "Genus")
  
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
                          number = 10, repeats=3,    # 3 repeats of 10-fold cross-validation, it will perform 10-fold cross-validation on the training data 3 times
                          summaryFunction=prSummary,
                          classProbs=T,
                          savePredictions = T,
                          verboseIter = TRUE)
    
    
    # Training ML Genus and Smoking habits
    # Running GBM here
    gbm <- train(Sex ~., 
                 data = df, 
                 method = "gbm",
                 metric = "AUC",
                 bag.fraction = 0.7, # bag.fraction (Subsampling fraction) - the fraction of the training set observations randomly selected to propose the next tree in the expansion. In this case, it adopts stochastic gradient boosting strategy. By default, it is 0.5. That is half of the training sample at each iteration. You can use fraction greater than 0.5 if training sample is small.
                 trControl = cntrl)
    
    # Running Random Forest
    rf <- train(Sex ~., 
                data = df, 
                method = "rf",
                metric = "AUC",
                trControl = cntrl)
    
    
    
    ## Analyze the ML models with  MLeval
    res <- evalm(list(gbm, rf), showplots = TRUE, title = paste0("Healthy - ", s),
                 gnames = c("GBM", "RF"),
                 cols = c("#ca0020", "#999999"))
    
    
    sc.list[[paste0("Healthy - ", s)]] <- res
    
  } 
 
  
  # My theme start here
  my_theme <- theme(legend.position = c(.65,.15),
                    legend.text = element_text(size=9),
                    legend.key = element_rect(fill = NA),
                    legend.background=element_blank(),
                    text = element_text(family = "Avenir"))

   # Periodontitis - Saliva
   x <- cbind(table(sc.list[["Periodontitis - Saliva"]][["proc"]][["data"]][["obs"]]))
   x <- as.character(paste(rownames(x), x, collapse = '\n'))
   x <- gsub('female', 'F:', x)
   x <- gsub('male', 'M:', x)
  p1 <- sc.list$`Periodontitis - Saliva`$roc +
          annotate("text", label = x, x=0, y=.85, vjust = 0, hjust=0, size = 3.5)+
          my_theme 
  # Set my color here
  p1 <- p1 + scale_color_manual(values = c("#b3494e", "#ff8288"))
  
  # Periodontitis - Plaque
  x <- cbind(table(sc.list[["Periodontitis - Plaque"]][["proc"]][["data"]][["obs"]]))
  x <- as.character(paste(rownames(x), x, collapse = '\n'))
  x <- gsub('female', 'F:', x)
  x <- gsub('male', 'M:', x)
  p2 <- sc.list$`Periodontitis - Plaque`$roc +
            annotate("text", label = x, x=0, y=.85, vjust = 0, hjust=0, size = 3.5)+
             my_theme 
  # Set my color here
  p2 <- p2 + scale_color_manual(values = c("#b3494e", "#ff8288"))
  
  # Periodontitis - Subgingival
  x <- cbind(table(sc.list[["Periodontitis - Subgingival"]][["proc"]][["data"]][["obs"]]))
  x <- as.character(paste(rownames(x), x, collapse = '\n'))
  x <- gsub('female', 'F:', x)
  x <- gsub('male', 'M:', x)
  p3 <- sc.list$`Periodontitis - Subgingival`$roc +
        annotate("text", label = x, x=0, y=.85, vjust = 0, hjust=0, size = 3.5)+
         my_theme 
  # Set my color here
  p3 <- p3 + scale_color_manual(values = c("#b3494e", "#ff8288"))
  
  
  
  
  # Healthy - Saliva
  x <- cbind(table(sc.list[["Healthy - Saliva"]][["proc"]][["data"]][["obs"]]))
  x <- as.character(paste(rownames(x), x, collapse = '\n'))
  x <- gsub('female', 'F:', x)
  x <- gsub('male', 'M:', x)
  p4 <- sc.list$`Healthy - Saliva`$roc +
    annotate("text", label = x, x=0, y=.85, vjust = 0, hjust=0, size = 3.5)+
    my_theme 
  # Set my color here
  p4 <- p4 + scale_color_manual(values = c("#017332", "#04bf55"))
  
  
  # Healthy - Plaque
  x <- cbind(table(sc.list[["Healthy - Plaque"]][["proc"]][["data"]][["obs"]]))
  x <- as.character(paste(rownames(x), x, collapse = '\n'))
  x <- gsub('female', 'F:', x)
  x <- gsub('male', 'M:', x)
  p5 <- sc.list$`Healthy - Plaque`$roc +
    annotate("text", label = x, x=0, y=.85, vjust = 0, hjust=0, size = 3.5)+
    my_theme 
  # Set my color here
  p5 <- p5 + scale_color_manual(values = c("#017332", "#04bf55"))
  
  # Healthy - Subgingival
  x <- cbind(table(sc.list[["Healthy - Subgingival"]][["proc"]][["data"]][["obs"]]))
  x <- as.character(paste(rownames(x), x, collapse = '\n'))
  x <- gsub('female', 'F:', x)
  x <- gsub('male', 'M:', x)
  p6 <- sc.list$`Healthy - Subgingival`$roc +
    annotate("text", label = x, x=0, y=.85, vjust = 0, hjust=0, size = 3.5)+
    my_theme 
  # Set my color here
  p6 <- p6 + scale_color_manual(values = c("#017332", "#04bf55"))
  


library(patchwork)
AUCROC <-  (p4 + p5 +p6) / (p1 + p2 + p3)
AUCROC
  
  ggsave(plot = AUCROC, 
         device = cairo_pdf,
         width = 8, 
         height = 11, 
         units = "in",
         filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/AUROC Figure_2E.pdf")
  
  
  

  
  
  
  
  
  
  
  
    