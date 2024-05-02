library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(ANCOMBC)
library(dplyr)

# Load PSM healthy dataset
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# Agglomerate taxa of the same type
pseq_phylum <- phyloseq::tax_glom(psm_healthy, taxrank = "Phylum")

# Define my vars for for cycle
MySites <- unique(pseq_phylum@sam_data$Site)

# Prepare a datafrafe for storing resuts
results <- data.frame()

# Start cycle here:
  for (i in MySites) {     
print(i)
      # Subsetting according to my strata
      foo <- subset_samples(pseq_phylum, Site == i)

      print(foo)
      cat("\n")
      print(table(foo@sam_data$Site, foo@sam_data$Sex))
      
      # Running ANCOM-BC
      set.seed(101)
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
      
      # Put stratas into output
      out$res$Site <- list(Site)

      # Storing results
      res <- data.frame(out$res)
      names(res) <- c("beta", "se", "W", "p_val", "q_val", "diff_abn", "Site")
      
      # Report OTU name fron rownames to dedicated coloum
      res$OTUs <- rownames(res)
      
      # Appending results
      results <- rbind(results, res)
      
      # make some cleaning
     # rm(out, res, foo)
      
      # Print on the screen
      cat(paste0("Fineshed ", i,  " - ", "\n\n\n"))
    }
  

# Print Phyla that are different
DifferentPhyla <- subset(results, diff_abn == TRUE)
DifferentPhyla 

# Show taxa
tax_table(psm_healthy)[,"Genus"]

# Defining Taxa of interest
MyPhyla <- c("Patescibacteria","Desulfobacterota", "Chloroflexi", "Deinococcota")
MyClass <- c("Saccharimonadia", "Desulfobulbia", "Anaerolineae", "Deinococci")

# My dataset with phyla of interest
PhylaOfInterest <- subset_taxa(psm_healthy, Phylum %in% MyPhyla)

# Absolute abundances
Absolute = psmelt(PhylaOfInterest)

# Relative abundance
Relative = transform_sample_counts(PhylaOfInterest, function(x){x / sum(x)})
Relative = psmelt(Relative)


library(ggpubr)
ggboxplot(data = Relative, x = "Sex", y = "Abundance", fill = "Sex",
          add = "jitter",
          facet.by = c("Site", "Phylum"))+
  stat_compare_means()

library(plyr)
cdata <- ddply(Relative, c("Site", "Phylum", "Sex", "Smoking"), summarize,
               N = length(Run))
cdata

