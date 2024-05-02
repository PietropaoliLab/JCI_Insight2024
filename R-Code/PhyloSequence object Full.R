library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(phyloseq)
library(microbiome)

# Load OTU table and TAXA
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/V3 V4 V3V4 SeqTable and Taxa complete.RData")

# Load MetaData
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Metadata for analisys.RData")

# Capture sample ID from OTU table
MyRun <- rownames(seqtable_full)
MyRun
MyRun <- gsub(".fastq", "", MyRun)
rownames(seqtable_full) <- MyRun
MyRun
rownames(MetaData) <- MetaData$Run

ps <- phyloseq(otu_table(seqtable_full, taxa_are_rows=FALSE), 
               sample_data(MetaData), 
               tax_table(taxa))


summarize_phyloseq(ps)
# Show available ranks in the dataset
rank_names(ps)

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# The following ensures that features with ambiguous phylum annotation are also removed. 
# Note the flexibility in defining strings that should be considered ambiguous annotation.
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))


tot <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
tot <- subset(tot, `2` > 5)
tot


# Define phyla to filter
filterPhyla = tot$Phylum
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, Phylum %in% filterPhyla)
ps1
ps



# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.03, alpha = 0.5, linetype = 2, color = "red") +  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


# Define prevalence threshold as 3% of total samples
prevalenceThreshold = 0.03 * nsamples(ps)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
ps2

# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))

ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
ps3

ps2
ps3


plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes", "Bacteroidota"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Sex",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}


# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})

plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)

plotBefore
plotAfter



psOrd = subset_taxa(ps3ra, Phylum == "Bacteroidota")
plot_abundance(psOrd, Facet = "Order", Color = NULL)+geom_boxplot(color = "green")

ps3
# Keep samples with non-zero total taxa
pruned <- prune_samples(sample_sums(ps3) > 0, ps3)
pruned

ps <- pruned


tax <- data.frame(tax_table(ps))

library(stringr)
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean$kiss <- NA

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fille holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  
  ## Fill in missing taxonomy
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax.clean <- tax.clean[,c(1:6,8)]               


n_seqs <- seq(1, nrow(tax.clean))
len_n_seqs <- nchar(max(n_seqs))

taxa_names <- paste("OTU", formatC(n_seqs, 
                                   width = len_n_seqs, 
                                   flag = "0"), sep = "_")

taxa_names <- paste0(taxa_names, "_", tax.clean$Species)

rownames(tax.clean) <- taxa_names
tax.clean <- tax.clean[1:6]

tax.clean <- as.matrix(tax.clean)

taxa_names(ps) <- taxa_names

head(taxa_names(ps))




save(ps, file = "Microbiome/Gender features of periodontal microbiome/Data/PlyloSeq Object for analysys.RData")


