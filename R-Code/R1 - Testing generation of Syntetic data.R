library(phyloseq)
library(here)
library(dplyr)
library(tidyr)
library(synthpop)


# Load Phyloseq
load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS_for_revision.RData"))
table(wgs@sam_data$Disease, wgs@sam_data$Eth)

# To simplify the WGS data
ps = (wgs)

# Select only "Caucasian" individuals
ps <- subset_samples(ps, Eth == "Caucasian")


# Check Disease vs Eth
table(ps@sam_data$Disease, ps@sam_data$Eth)

total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_genus = transform_sample_counts(ps, standf)


# Sample Data has Sample ID per row
normalized_genus = tax_glom(normalized_genus, taxrank = "Genus", NArm = TRUE)
normalized_genus


normalized_genus_RelAbn = transform_sample_counts(normalized_genus, function(OTU) OTU/sum(OTU) * 100)
normalized_genus_RelAbn


genus <- psmelt(normalized_genus_RelAbn)


# Keep anly Genus with Relative abuntance >= 5%
genus <- subset(genus, Abundance >= 5)


# Pivot the data frame
pivot_genus <- genus %>%
  pivot_wider(id_cols = c(Sample, Eth, Sex, Disease), names_from = OTU, values_from = Abundance, values_fill = 0)


# Keep anly Caucasians
pivot_genus <- subset(pivot_genus, Eth == "Caucasian")
pivot_genus <- data.frame(pivot_genus)

rownames(pivot_genus) <- paste(seq(1:nrow(pivot_genus)),pivot_genus$Sample, pivot_genus$BioSample, sep = "_")
pivot_genus <- pivot_genus[4:58]


codebook.syn(pivot_genus)

