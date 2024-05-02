library(phyloseq)
# Load my 16S data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Load NHANES III
load("~/Google Drive/01- UNIVERSITA/AQ Projects/NHANES/NHANES III/data/NHANES for analysis DEC2018.RData")

# Working on microbiome data
# PhyloSeq metadata
MetaDataPerio <- data.frame(psm_perio@sam_data)
MetaData <- MetaDataPerio[, c("Run", "Sex", "AgeYrs", "Periodontitis", "Site", "Smoking")]

# Subsetting only subgingival
MetaData <- subset(MetaData, Site == "Subgingival")

# After subsetting unknown smoking habits
MetaData <- subset(MetaData, !is.na(Smoking))

MetaData <- MetaData[-5]

# Screening data
table(MetaData$Periodontitis, MetaData$Sex)


# Working on NHANES III
nhanes <- nhanes[,c("SEQN", "SEX", "AGE",  "periodontitis", "AT.LEAST.100.CIGS")]

# Subsetting only individuals with valis periodontal classification
# starting from 8153
# Valid periodontal clasification: 5825
nhanes <- subset(nhanes, !is.na(periodontitis))

# SMOKINNG
nhanes <- subset(nhanes, !is.na(AT.LEAST.100.CIGS))

# Changing colnames for consistency
names(nhanes) <- c("SEQN", "Sex", "AgeYrs", "Periodontitis", "Smoking")

# Verifining smoking var
table(nhanes$Smoking)

# Same factors for consistency
unique(nhanes$Sex)
nhanes$Sex <- as.character(nhanes$Sex)
nhanes$Sex <- gsub("Male", "male", nhanes$Sex)
nhanes$Sex <- gsub("Female", "female", nhanes$Sex)
nhanes$Periodontitis <- ifelse(nhanes$Periodontitis == 1, "Periodontitis", "Healthy")
nhanes$Smoking <- ifelse(nhanes$Smoking == 1, "Yes", "No")


# Matching pairs
library(tidyverse)

MetaData$NHANES_Pairs <- NA

for (i in 1:nrow(MetaData)) {
  
  small.attributes <- select(MetaData[i, ], -Run, -NHANES_Pairs)
  
  big_matches <- nhanes %>% 
    mutate(across(-SEQN, ~.x == small.attributes[[cur_column()]])) %>% 
    filter(if_all(-SEQN, ~. == T) & !(SEQN %in% MetaData$NHANES_Pairs)) 
  MetaData$NHANES_Pairs[i] <- big_matches$SEQN[1]
  
}
rm(big_matches, i, small.attributes)

Matched <- filter(nhanes, SEQN %in% MetaData$NHANES_Pairs)

names(MetaData)[6] <- "SEQN"

Pairs <- left_join(MetaData, nhanes, by ="SEQN", suffix = c("_16S", "_NHANES"))
Pairs <- na.omit(Pairs)

# Subsetting dataframes accoring to pairs
Perio <- subset(Pairs, Periodontitis_16S == "Periodontitis")

# Periodontitis pairs
MyRun <- Perio$Run
PerioMicrobiome <- subset_samples(psm_perio, Run %in% MyRun)



# NHANES pairs
# Load NHANES III
load("~/Google Drive/01- UNIVERSITA/AQ Projects/NHANES/NHANES III/data/NHANES for analysis.RData")
nhanes <- subset(nhanes, SEQN %in% Pairs$SEQN)
nhanes <- nhanes[,c("SEQN", "SEX", "AGE", "periodontitis","P.gingivalis", 
                    "P.intermedia", "P.nigrescens", "T.forsythia", "A.actino-mix", 
                    "F.nucleatum", "S.oralis",  "M.micros", "C.rectus", "E.corrodens", "E.nodatum", "S.intermedius", 
                    "C.ochracea", "V.parvula", "A.naeslundii", "P.melaninogenica", 
                    "S.noxia", "T.denticola", "S. mutans")]
names(nhanes) <- gsub(" ", "", names(nhanes))
names(nhanes) <- gsub("-", "_", names(nhanes))


nhanes$SEX <- as.character(nhanes$SEX)
nhanes$SEX <- gsub("Male", "male", nhanes$SEX)
nhanes$SEX <- gsub("Female", "female", nhanes$SEX)
nhanes$periodontitis <- ifelse(nhanes$periodontitis == 1, "Periodontitis", "Healthy")

# Performing Log10 normalization
MyStrains <- names(nhanes[5:23])
Normalized <- data.frame("SEQN" = nhanes$SEQN)

for (i in MyStrains) {
  foo <- log10(nhanes[i]+1) # add +1 to avoid -Inf results
  names(foo) <- paste0(names(foo),".log10")
  Normalized <- cbind(Normalized, foo)
  rm(foo)
}
rm(MyStrains)

# Perform RowSums for relative abundance
Normalized$AbSums <- rowSums(Normalized[2:20])

# Performing relative abundance
RelAbn <- data.frame("SEQN" = Normalized$SEQN)
MyStrains <- names(Normalized[2:20])

for (i in MyStrains) {
  foo <- Normalized[i]/Normalized$AbSums
  names(foo) <- gsub(".log10", ".RelAbn", names(foo))
  RelAbn <- cbind(RelAbn, foo)
  rm(foo)
}

# Merging dataframes
nhanes <- merge(nhanes, Normalized, by = "SEQN")
nhanes <- merge(nhanes, RelAbn, by = "SEQN")

## Generating Genus Ab in NHANES
# Eubacterium = E.nodatum
nhanes$`[Eubacterium] nodatum group` <- nhanes$E.nodatum.RelAbn

# Actinomyces = A.naeslundii
nhanes$Actinomyces <- nhanes$A.naeslundii.RelAbn

# Aggregatibacter = a.actimo.mix
nhanes$Aggregatibacter <- nhanes$A.actino_mix.RelAbn

# Campylobacter = C.rectus
nhanes$Campylobacter <- nhanes$C.rectus.RelAbn

# Capnocytophaga = C.ocracea
nhanes$Capnocytophaga <- nhanes$C.ochracea.RelAbn

# Eikenella = E.corrodens
nhanes$Eikenella <- nhanes$E.corrodens.RelAbn

# Fusobacterium = F.nucleatum
nhanes$Fusobacterium <- nhanes$F.nucleatum.RelAbn

# Parvimonas = M.micros
nhanes$Parvimonas <- nhanes$M.micros.RelAbn

# Porphyromonas = P.gingivalis
nhanes$Porphyromonas <- nhanes$P.gingivalis.RelAbn

# Prevotella = P.intermedia; P.nigrescens; P.melaninogenica
nhanes$Prevotella <- rowSums(nhanes[,c("P.intermedia.RelAbn", "P.nigrescens.RelAbn")])

# Prevotella_7 = S.noxia
nhanes$Prevotella_7 <- nhanes$P.melaninogenica.RelAbn

# Selenomonas = S.noxia
nhanes$Selenomonas <- nhanes$S.noxia.RelAbn

# Streptococcus = S.oralis; S.intermedius; S.mutans
nhanes$Streptococcus <- rowSums(nhanes[,c("S.oralis.RelAbn", "S.intermedius.RelAbn", "S.mutans.RelAbn")])

# Tannerella = T.forsythia
nhanes$Tannerella <- nhanes$T.forsythia.RelAbn

# Treponema = T.denticola
nhanes$Treponema <- nhanes$T.denticola.RelAbn

# Veillonella = V.parvula
nhanes$Veillonella <- nhanes$V.parvula.RelAbn


## Generating Phylum Ab in NHANES
# Actinobacteriota = Actinomyces
nhanes$Actinobacteriota <- nhanes$Actinomyces

# Bacteroidetes = Porphyromonas; Prevotella; Tannerella; 
nhanes$Bacteroidetes <- rowSums(nhanes[,c("Porphyromonas", "Prevotella", "Tannerella")])

# Bacteroidota = Capnocytophaga; Prevotella_7; 
nhanes$Bacteroidota <- rowSums(nhanes[,c("Capnocytophaga", "Prevotella_7")])

# Campylobacterota = Campylobacter
nhanes$Campylobacterota <- nhanes$Campylobacter

# Firmicutes = Streptococcus; Parvimonas; `[Eubacterium] nodatum group`; Veillonella; Selenomonas; 
nhanes$Firmicutes <- rowSums(nhanes[,c("Parvimonas", "[Eubacterium] nodatum group", "Veillonella", "Selenomonas")])

# Fusobacteriota = Fusobacterium
nhanes$Fusobacteriota <- nhanes$Fusobacterium

# Proteobacteria = Aggregatibacter; Eikenella; 
nhanes$Proteobacteria <- rowSums(nhanes[,c("Aggregatibacter",  "Eikenella")])

# Spirochaetota = Treponema
nhanes$Spirochaetota <- nhanes$Treponema


# Some cleaning
rm(Normalized, RelAbn, MyRun, i, Matched, MetaData, Perio, MyStrains)



# Subsetting TAXA according to Phylum in NHANES antibodies
MyPhylum <- c("Bacteroidetes", "Proteobacteria", "Fusobacteriota", "Firmicutes", 
              "Campylobacterota", "Bacteroidota", "Actinobacteriota", "Spirochaetota")

MyGenus <- c("Porphyromonas", "Prevotella", "Tannerella", "Aggregatibacter", 
             "Fusobacterium", "Streptococcus", "Parvimonas", "Campylobacter", 
             "Eikenella", "[Eubacterium] nodatum group", "Capnocytophaga", 
             "Veillonella", "Actinomyces", "Prevotella_7", "Selenomonas", 
             "Treponema")


# Selecting only TAX of interest
PhylumNHANES <- nhanes[,c("SEQN", MyPhylum)]
PhylumNHANES$Tax <-"Phylum"

GenusNHANES <-  nhanes[,c("SEQN", MyGenus)]
GenusNHANES$Tax <-"Genus"

x <- Pairs[,c("Run", "SEQN")]

GenusNHANES <- merge(GenusNHANES, x)
PhylumNHANES <- merge(PhylumNHANES, x)

# My Phyla antibodies to be sum
# MyAb <- c("Bacteroidota", "Firmicutes", "Fusobacteriota", "Proteobacteria", "Spirochaetota", "Actinobacteriota") # Phylum
MyAb <- c("Tannerella","Streptococcus","Fusobacterium","Treponema","Campylobacter","Prevotella","Prevotella_7") # Genus


# Sum my ab of interest
PhylumNHANES$AbTitersOfInterest <- rowSums(PhylumNHANES[MyAb], na.rm = TRUE)
PhylumNHANES <- PhylumNHANES[,c("SEQN", "Run","AbTitersOfInterest")]

GenusNHANES$AbTitersOfInterest <- rowSums(GenusNHANES[MyAb], na.rm = TRUE)
GenusNHANES <- GenusNHANES[,c("SEQN", "Run","AbTitersOfInterest")]


####################################
#
##  Microgenderome Index - MGI
#
####################################
tmp <- subset_samples(psm_perio, Run %in% PhylumNHANES$Run)
tmp <- tax_glom(tmp, taxrank = "Genus")
tmp <- microbiome::transform(tmp, transform = "log10p")
tmp <- microbiome::transform(tmp, transform = "compositional")


# Subset Genera of interest
MaleGenus <- c("Olsenella", "Peptoanaerobacter")
FemaleGenus <- c("Streptococcus", "Fusobacterium", "Tannerella", "Treponema", 
                 "Campylobacter", "Prevotella")
Genus <- c(MaleGenus, FemaleGenus)


# Prepare data for MGI
df <- subset_taxa(tmp, Genus %in% Genus)
df <- data.frame(otu_table(df))
colnames(df) <- gsub(".*Genus_", "", colnames(df))

# Sums relative abundance of enriched genus in male or in female
df$MaleEnriched <- rowSums(df[MaleGenus])
df$FemaleEnriched <- rowSums(df[FemaleGenus])

# Compute MGI
df$MGI <- log((1+df$FemaleEnriched)/(1+df$MaleEnriched))

# Store MGI into tmp phylobject
tmp@sam_data$MGI <- df$MGI

tmp <- subset_samples(tmp, Run %in% GenusNHANES$Run)
tmp <- sample_data(tmp)
tmp <- tmp[,c("Run","Sex", "Periodontitis", "Smoking", "MGI")]
tmp <- data.frame(tmp)

MyFinalDataset <- merge(GenusNHANES, tmp, by= "Run")


rm("df", "FemaleGenus", "Genus", "GenusNHANES", "MaleGenus", "MetaDataPerio", 
    "MyAb", "MyGenus", "MyPhylum", "nhanes", "Pairs", 
    "PerioMicrobiome", "PhylumNHANES", "psm_perio", "tmp", "x")

library(ggpubr)
ggscatter(subset(MyFinalDataset, Smoking == "No"), 
          x = "MGI", y =  "AbTitersOfInterest",
          color = "Sex",
          add = "reg.line",                                       # Add regression line
          conf.int = TRUE,                                        # Add confidence interval
          fullrange = TRUE,
          shape = 21,
          size = 3)+
  stat_cor(aes(color = Sex),method = "pearson", label.x = c(0.1, 0.1), label.y = c(0.44, 0.45))


