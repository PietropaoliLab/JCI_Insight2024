library(phyloseq)
library(synthpop)
library(dplyr)
library(tidyr)
library(here)

load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData"))
load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData"))

# Define my subseto of interest
perio <- subset_samples(psm_perio, Site == "Subgingival" & !is.na(Smoking))
healty <- subset_samples(psm_healthy, Site == "Subgingival" & !is.na(Smoking))

# Transform in compositional
# perio <- microbiome::transform(perio, transform = "compositional")
# healty <- microbiome::transform(healty, transform = "compositional")

# melting
perio <- psmelt(perio)
healty <- psmelt(healty)


# perio$OTU <- gsub("^OTU_[0-9]+_", "", perio$OTU)
# healty$OTU <- gsub("^OTU_[0-9]+_", "", healty$OTU)


# Pivot the data frame
pivot_perio <- perio %>%
                        pivot_wider(id_cols = c(Sample, Sex, Perio, Smoking, AgeYrs), names_from = OTU, values_from = Abundance, values_fill = 0)
# Pivot the data frame
pivot_healty <- healty %>%
                        pivot_wider(id_cols = c(Sample, Sex, Perio, Smoking, AgeYrs), names_from = OTU, values_from = Abundance, values_fill = 0)



codebook.syn(pivot_perio[-1])
codebook.syn(pivot_healty[-1])


# Note which variables have missing values, especially those that are not coded as the R missing value NA.  
# For example the value -9 often signifies missing data for positive items like income. 
# These can be identified to the syn() function via the cont.na parameter.
# Note any variables that ought to be derivable from others, e.g. discharge date from length-of-stay and 
# admission date. These could be omitted and recalculated after synthesis or calculated as part of the synthesis 
# process by setting their method to passive (see ?syn.passive).
# Also note any variables that should obey rules that depend on other variables. 
# For example, the number of cigarettes smoked should be zero or missing for non-smokers. 
# You can set the rules with the parameters rules  and rvalues of the syn() function. 
# The syn() function will warn you if the rule is not obeyed in the observed data.

# You are now ready to do your first synthesis
PerioSyn <- syn(pivot_perio[-1], seed = 101, k = 100) # a size of the synthetic data set (k x p), which can be smaller or greater than the size of the original data set (n x p). The default is nrow(data) which means that the number of individuals in the synthesised data is the same as in the original (observed) data (k = n).
HealtySyn <- syn(pivot_healty[-1], seed = 101, k = 100) # a size of the synthetic data set (k x p), which can be smaller or greater than the size of the original data set (n x p). The default is nrow(data) which means that the number of individuals in the synthesised data is the same as in the original (observed) data (k = n).


# To get an overview, use the summary() function for a synds object, e.g.
summary(PerioSyn)
summary(HealtySyn)


# To do an initial comparison of the original and synthetic data as tables and histograms use the compare() function, e.g.
# compare(PerioSyn, pivot_perio[-1], stat = "percents")
# compare(HealtySyn, pivot_healty[-1], stat = "percents")


PerioSyn <- PerioSyn$syn
HealtySyn <- HealtySyn$syn

Syn <- rbind(PerioSyn, HealtySyn)

table(Syn$Perio, Syn$Sex)

tax_table(psm_healthy)

# taxa table are identical?
identical(tax_table(psm_perio), tax_table(psm_healthy))

# producing objects for Syntetic Phyloseq
# Taxonomy table
TAX <- tax_table(psm_perio)

# OTU table
OTU <- as.matrix(t(Syn[4:137]))
OTU = otu_table(OTU, taxa_are_rows = TRUE)

# metadata
Meta <- sample_data(Syn[1:4])


# Syntetic Phyloseq object creation
ps <- phyloseq(OTU, TAX, Meta)


# Count disease vs sex
table(ps@sam_data$Perio, ps@sam_data$Sex)


save(ps, file = here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Syntetic_16S_dataset.RData"))

rm(list = ls())
gc()