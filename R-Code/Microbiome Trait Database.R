# Exploring the functional composition of the human microbiome using 
# a hand-curated microbial trait database
# doi.org/10.1186/s12859-021-04216-2


library(readr)
trait <- read_csv("Microbiome/Gender features of periodontal microbiome/Trait Database.csv")

save(trait, file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Microbiome Trait.Rdata")
