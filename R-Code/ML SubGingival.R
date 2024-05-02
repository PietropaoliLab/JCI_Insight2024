library(phyloseq) 
library(microbiome)
library(caret)
library(MLeval)
library(reshape2)

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
healthy <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")
perio <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

healthy <- phyloseq::subset_samples(healthy, Site == "Subgingival")
perio <- phyloseq::subset_samples(perio, Site == "Subgingival")



# Melting Phyloseq Objects
healthy <- psmelt(healthy)
perio <- psmelt(perio)

# My vars for training ML
myvars<-c("Run","Sex", "Genus", "Abundance") 


# Define the DF of interest providing NA cleaning
df<-healthy[myvars]
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
             trControl = cntrl)

# Running Random Forest
rf <- train(Sex ~., 
            data = df, 
            method = "rf",
            metric = "AUC",
            trControl = cntrl)



## Analyze the ML models with  MLeval
res <- evalm(list(gbm, rf), showplots = TRUE, title = "Healthy",
             gnames = c("GBM", "RF"),
             cols = c("#ca0020", "#999999"))

evalm(gbm)
evalm(rf)

# My Prediction
HealthyPrediction <- res$probs


## get ROC
res$roc

## AUC and other metrics
res$optres


## VAriable importance into df
library(gbm)
VarImportanceGBM <- varImp(gbm, scale = TRUE)
VarImportanceRF <- varImp(rf, scale = TRUE)


V1 <- data.frame(VarImportanceGBM$importance, "Model" = "GBM")
V1$Genus <- rownames(V1)
V1 <- V1[order(-V1$Overall), ]


V2 <- data.frame(VarImportanceRF$importance, "Model" = "RF")
V2$Genus <- rownames(V2)
V2 <- V2[order(-V2$Overall), ]


VarImportanceHealthy <- rbind(V1[1:20,], V2[1:20,])


library(ggplot2)
VarImportanceHealthy <- VarImportanceHealthy[order(-VarImportanceHealthy$Overall), ]


MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = c(0.85, 0.1),
        legend.key = element_blank(),
        aspect.ratio = 2/1,
        legend.background=element_blank(),
        plot.margin=unit(c(0.5,0,0.5,0), "cm"))

p <- ggplot(VarImportanceHealthy, aes(x= reorder(Genus, Overall), y = Overall, group = Model, fill = Model))+
  geom_col(position = position_dodge2(preserve = "single"), width = .5) +
  scale_fill_manual(values = c("#ca0020", "#999999"))+
  labs(x="", y="Feature importance (%)", title = "Healthy")+
  coord_flip() +
  MyTheme

############################################################################################



# Define the DF of interest providing NA cleaning
df<-perio[myvars]
df<-na.omit(df)

# Generate abundance table for ML
df$AbundanceLog10 <- log10(df$Abundance + 1)

abn <- df[c("Run", "Genus", "AbundanceLog10")]
abn$Genus <- gsub(" ", "_", abn$Genus)
abn$Genus <- gsub("-", "_", abn$Genus)
abn$Genus <- gsub("\\[", "", abn$Genus)
abn$Genus <- gsub("\\]", "", abn$Genus)
abn$Genus <- gsub("\\.", "", abn$Genus)

abn_wide <- dcast(abn, Run ~ Genus, value.var="AbundanceLog10")

df <- df[c("Run", "Sex")]
df <- unique(df)

df <- merge(df, abn_wide, by = "Run")
rownames(df) <- df$Run

df$Sex <- factor(df$Sex)
#df$Smoking <- factor(df$Smoking)
#df$Site <- factor(df$Site)


df <- na.omit(df)

df <- df[-1]

## partition for test and train at 80% train
set.seed(101)
cntrl <- trainControl(method = "repeatedcv",
                      number = 10, repeats = 3,
                      summaryFunction=prSummary,
                      classProbs=T,
                      savePredictions = T,
                      verboseIter = TRUE)


# Training ML Genus and Smoking habits
gbm <- train(Sex ~., 
             data = df, 
             method = "gbm",
             metric = "AUC",
             trControl = cntrl)
rf <- train(Sex ~., 
            data = df, 
            method = "rf",
            metric = "AUC",
            trControl = cntrl,
            verbose = TRUE)



## run MLeval

res1 <- evalm(list(gbm, rf), showplots = TRUE, title = "Periodontitis",
              gnames = c("GBM", "RF"),
              cols = c("#ca0020", "#999999"))

evalm(gbm)
evalm(rf)

# My Periodontitis Predictions
PerioPrediction <- res1$prob

## get ROC
res1$roc

## AUC and other metrics
res1$optres


## VAriable importance into df
library(gbm)
VarImportanceGBM <- varImp(gbm, scale = TRUE)
VarImportanceRF <- varImp(rf, scale = TRUE)


V1 <- data.frame(VarImportanceGBM$importance, "Model" = "GBM")
V1$Genus <- rownames(V1)
V1 <- V1[order(-V1$Overall), ]


V2 <- data.frame(VarImportanceRF$importance, "Model" = "RF")
V2$Genus <- rownames(V2)
V2 <- V2[order(-V2$Overall), ]


VarImportancePerio <- rbind(V1[1:20,], V2[1:20,])


library(ggplot2)
VarImportancePerio <- VarImportancePerio[order(-VarImportancePerio$Overall), ]


MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        aspect.ratio = 2/1,
        legend.background=element_blank(),
        plot.margin=unit(c(0.5,0,0.5,0), "cm"))

p1 <- ggplot(VarImportancePerio, aes(x= reorder(Genus, Overall), y = Overall, group = Model, fill = Model))+
  geom_col(position = position_dodge2(preserve = "single"), width = .5) +
  scale_fill_manual(values = c("#ca0020", "#999999"))+
  labs(x="", y="Feature importance (%)", title = "Periodontitis")+
  coord_flip() +
  MyTheme+
  theme(plot.margin=unit(c(0.5, 0.5, 0.5,-0.5), "cm"),
        legend.position = "none")
p1

library(gridExtra)
library(grid)
grid.draw(cbind(ggplotGrob(p), ggplotGrob(p1), size="last"))

