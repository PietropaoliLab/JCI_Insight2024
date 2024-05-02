library(microbiomeSeq)  #load the package
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(showtext)
library(patchwork)
library(here)
library(readr)
library(rentrez)
library(XML)
library(readxl)
library(tidyverse)
loadfonts(device = "pdf")

# Load Phyloseq
load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS_for_revision.RData"))

# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 1/3,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

# To simplify the WGS data
ps = (wgs)

total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_ps = transform_sample_counts(ps, standf)


# Sample Data has Sample ID per row
normalized_ps_phylum1 = tax_glom(normalized_ps, "Phylum", NArm = TRUE)
normalized_ps_phylum1


normalized_ps_phylum1_relabun = transform_sample_counts(normalized_ps_phylum1, function(OTU) OTU/sum(OTU) * 100)


df <- psmelt(normalized_ps_phylum1_relabun)
ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Relative Abundance", title = "Phylum Relative Abundance")+
  facet_wrap(Sex ~ Disease, scales = "free", drop = TRUE) +
  MyTheme

names(df)

x <- df[,c("Sample","Site","Geo", "Eth", "Sex", "Disease","Phylum","Abundance" )]
# Removing NAs 
x <- na.omit(x)
length(unique(x$Sample))
table(x$Sex, x$Disease)

library(dplyr)
library(tidyr)

# Pivot the data frame
pivot_df <- x %>%
  pivot_wider(names_from = Phylum, values_from = Abundance, values_fill = 0)


pivot_df <- subset(pivot_df, Eth == "Caucasian")
pivot_df <- pivot_df[,c(1, 5:20)]
pivot_df <- pivot_df[,c(2:13,15)]

write_csv(pivot_df, "~/Downloads/For_shyntetic_data_generation.csv")

library(synthpop)
codebook.syn(pivot_df)

# Use the output to do the following things to make your data ready to be synthesised:
# Remove any identifiers, e.g. study number.
# Change any character (text) variables into factors and rerun codebook.syn() after this. 
# The syn() function will do this conversion for you but it is better that you do it first.
pivot_df$Sex <- factor(pivot_df$Sex)
pivot_df$Disease <- factor(pivot_df$Disease)


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
mysyn <- syn(pivot_df, seed = 101, k = 1000) # a size of the synthetic data set (k x p), which can be smaller or greater than the size of the original data set (n x p). The default is nrow(data) which means that the number of individuals in the synthesised data is the same as in the original (observed) data (k = n).

# To get an overview, use the summary() function for a synds object, e.g.
summary(mysyn)

# To do an initial comparison of the original and synthetic data as tables and histograms use the compare() function, e.g.
compare(mysyn, pivot_df, stat = "percents")

sys <- mysyn$syn



