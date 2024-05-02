library(phyloseq)
library(tableone)
library(here)

# Load Phyloseq
load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS_for_revision.RData"))

# My MetaData variales
MyVars <- c( "LibraryStrategy", "LibrarySelection", 
            "LibrarySource", "LibraryLayout", 
            "Platform", "Model", 
            "CenterName", "Site", "Geo", "Eth", "Sex", "Age", "Disease", 
            "Smoke", "BioProject")

df <-  data.frame(wgs@sam_data)
tableOne <- CreateTableOne(vars = MyVars, strata = "Disease", data = df)

print(tableOne, quote = TRUE, noSpaces = TRUE, showAllLevels = TRUE, missing = TRUE)
