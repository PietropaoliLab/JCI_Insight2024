library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ANCOMBC)
library(dplyr)

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM PlyloSeq Object for analysys.RData")

pseq_genus <- phyloseq::tax_glom(psm_ps, taxrank = "Genus")

# ####################################

out = ancombc(phyloseq = pseq_genus, 
              formula = "Sex + Site + Periodontitis + Smoking", 
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

# ####################################

res = out$res

# Beta results
tab_coef = res$beta
colnames(tab_coef) = c("Beta_Sex", "Beta_Saliva", "Beta_Subgingival", "Beta_Periodontitis", "Beta_Smoking")
tab_coef$OTU <- rownames(tab_coef)

# SE results
tab_se = res$se
colnames(tab_se) = c("SE_Sex", "SE_Saliva", "SE_Subgingival", "SE_Periodontitis", "SE_Smoking")
tab_se$OTU <- rownames(tab_se)

# P valus results
tab_p = res$p_val
colnames(tab_p) = c("Pval_Sex", "Pval_Saliva", "Pval_Subgingival", "Pval_Periodontitis", "Pval_Smoking")
tab_p$OTU <- rownames(tab_p)

# P valus results
tab_q = res$q
colnames(tab_q) = c("Padj_Sex", "Padj_Saliva", "Padj_Subgingival", "Padj_Periodontitis", "Padj_Smoking")
tab_q$OTU <- rownames(tab_q)

tab_diff = res$diff_abn
colnames(tab_diff) = c("DiffAbn_Sex", "DiffAbn_Saliva", "DiffAbn_Subgingival", "DiffAbn_Periodontitis", "DiffAbn_Smoking")
tab_diff$OTU <- rownames(tab_diff)

ToMerge <- list(tab_coef, tab_se, tab_p, tab_q, tab_diff)

results <- Reduce(merge, ToMerge)

rm(tab_coef, tab_se, tab_p, tab_q, tab_diff, res, ToMerge)

MyTaxa <- results$OTU
df <- prune_taxa(MyTaxa, pseq_genus)
df <- as.data.frame(tax_table(df))
df$OTU <- rownames(df)

results <- left_join(results, df)

rm(df)












MyPerio <- unique(psm_ps@sam_data$Periodontitis)
MySites <- unique(psm_ps@sam_data$Site)
MySmoke <- c("No", "Yes")

results <- data.frame()

for (Perio in MyPerio) {
   for (Site in MySites) {
     for (Smoke in MySmoke) {
       
foo <- subset_samples(psm_ps, Periodontitis == MyPerio & Site == Site & Smoking == Smoke)


out = ancombc(phyloseq = foo, 
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

out$res$Site <- list(Site)
out$res$Perio <- list(Perio)
out$res$Smoking <- list(Smoke)


res <- data.frame(out$res)
names(res) <- c("beta", "se", "W", "p_val", "q_val", "diff_abn", "Site", "Periodontitis", "Smoking")

results <- rbind(results, res)

rm(res, foo)
cat(paste0("Fineshed ", Site, " - ", Perio, " - ", "Smoke=", Smoke, "\n"))
     }
   }
}


results$OTU <- rownames(results)
results <- subset(results, diff_abn == TRUE)
x <- as.character(results$OTU)
x <- gsub("OTU_", "", x)
x <- as.numeric(x)
x


y <- data.frame()
for (i in x) {
   print(i)
   foo <- as.data.frame(tax_table(pseq_genus)[i, 1:6])
   y <- rbind(y, foo)
   rm(foo)
   
}



df <- prune_taxa(x, psm_ps)

foo <- as.data.frame(tax_table(pseq_genus)[112, 1:6])

df <- prune_taxa(x, pseq_genus)
df <- as.data.frame(tax_table(df))

df$OTU <- rownames(df)

results <- left_join(results, df)

tax_table(pseq_genus)[tax_table(pseq_genus)[,"Order"]== "","Order"] <- "Unmatched_phylum"

