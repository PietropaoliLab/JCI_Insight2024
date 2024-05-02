
library("SIAMCAT")
library(phyloseq)
library(microbiome)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

data("feat_crc_zeller", package="SIAMCAT")
data("meta_crc_zeller", package="SIAMCAT")
head(meta.crc.zeller)

feat.crc.zeller <- data.frame(feat.crc.zeller) # OTU table - Rownames = OTUs and colnames = patienets ID
meta_crc_zeller <- data.frame(meta_crc_zeller) # Metedata  Rownames = patienets ID and colnames = Metadata vars names
head(meta_crc_zeller)

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")
pseq_genus <- microbiome::transform(pseq_genus, "compositional")

View(pseq_genus@otu_table)



meta <- pseq_genus@sam_data
meta <- data.frame(meta[,c("Sex", "AgeYrs", "Site", "BioProject", "size_MB")])
meta <- na.omit(meta)
meta$Sex <- factor(meta$Sex)

label.meta <- create.label(meta=meta,
                             label="Sex", case="female") 



sc.obj <- siamcat(feat=feat,
                  label=label.meta,
                  meta=meta)

show(sc.obj)

sc.obj <- filter.features(sc.obj,
                          filter.method = 'abundance',
                          cutoff = 0.001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)

 association.plot(sc.obj, sort.by = 'fc', 
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = '~/confounder_plots.pdf',
                   meta.in = NULL, feature.type = 'filtered')

sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                              norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

sc.obj <-  create.data.split(sc.obj, num.folds = 5, num.resample = 2)

sc.obj <- train.model(sc.obj, method = "randomForest")

sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
 
head(pred_matrix)
sc.obj <-  evaluate.predictions(sc.obj)
dev.off()
model.evaluation.plot(sc.obj)

model.interpretation.plot(sc.obj, fn.plot = '~/interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

