# Load metadating metadata
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/SampleMetadata_4VM.RData")

# Loading TAX
V3V4_tax <- readRDS("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Seq/V3V4/V3V4_tax.rds")

# Loading OTU
V3V4_seqtable <- readRDS("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Seq/V3V4/V3V4_seqtable.rds")


# Capture sample ID from OTU table
MyRun <- rownames(V3V4_seqtable) 

# Subsetting MetaData according to ID runs
df <- subset(SampleMetadata, Run %in% MyRun)

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

# Generating PyloSequence Object for furter Analysis
V3V4_ps <- phyloseq(otu_table(V3V4_seqtable, taxa_are_rows=FALSE), # OTU
                    sample_data(SampleMetadata), # Sample Metadata
                    tax_table(V3V4_tax))

V3V4_ps <- prune_samples(sample_names(V3V4_ps) != "Mock", V3V4_ps) # Remove mock sample

V3V4_ps <- ps_V3V4

plot_richness(V3V4_ps, x="Sex", measures=c("Shannon", "Simpson"), color="Perio")
