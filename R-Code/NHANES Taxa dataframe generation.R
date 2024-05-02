library(phyloseq)
# Load my 16S data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Load NHANES III
load("~/Google Drive/01- UNIVERSITA/AQ Projects/NHANES/NHANES III/data/NHANES for analysis.RData")

# Working on microbiome data
# PhyloSeq metadata
MetaDataHealthy <- data.frame(psm_healthy@sam_data)
MetaDataPerio <- data.frame(psm_perio@sam_data)
MetaDataHealthy <- MetaDataHealthy[, c("Run", "Sex", "AgeYrs", "Periodontitis", "Site")]
MetaDataPerio <- MetaDataPerio[, c("Run", "Sex", "AgeYrs", "Periodontitis", "Site")]
MetaData <- rbind(MetaDataHealthy, MetaDataPerio)

# Select only individuals with Subgingival microbiome
# Before subsetting 570 individuals
# After subsetting 199 
MetaData <- subset(MetaData, Site == "Subgingival")
MetaData <- MetaData[1:4]

# Screening data
table(MetaData$Periodontitis, MetaData$Sex)

# some cleaning
rm(MetaDataHealthy, MetaDataPerio)


# Working on NHANES III
nhanes <- nhanes[,c("SEQN", "SEX", "AGE",  "periodontitis")]

# Subsetting only individuals with valis periodontal classification
# starting from 8153
# Valid periodontal clasification: 5825
nhanes <- subset(nhanes, !is.na(periodontitis))

# Changing colnames for consistency
names(nhanes) <- c("SEQN", "Sex", "AgeYrs", "Periodontitis")

# Same factors for consistency
unique(nhanes$Sex)
nhanes$Sex <- as.character(nhanes$Sex)
nhanes$Sex <- gsub("Male", "male", nhanes$Sex)
nhanes$Sex <- gsub("Female", "female", nhanes$Sex)
nhanes$Periodontitis <- ifelse(nhanes$Periodontitis == 1, "Periodontitis", "Healthy")


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

names(MetaData)[5] <- "SEQN"

Pairs <- left_join(MetaData, nhanes, by ="SEQN", suffix = c("_16S", "_NHANES"))
Pairs <- na.omit(Pairs)

# Subsetting dataframes accoring to pairs
Perio <- subset(Pairs, Periodontitis_16S == "Periodontitis")
Healthy <- subset(Pairs, Periodontitis_16S == "Healthy")

# Periodontitis pairs
MyRun <- Perio$Run
PerioMicrobiome <- subset_samples(psm_perio, Run %in% MyRun)


# Healthy pairs
MyRun <- Healthy$Run
HealthyMicrobiome <- subset_samples(psm_healthy, Run %in% MyRun)


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
rm(Normalized, RelAbn, MyRun, i, Matched, MetaData, Perio, Healthy, MyStrains)



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

x <- Pairs[,c("Run", "SEQN","Sex_NHANES", "AgeYrs_NHANES", "Periodontitis_NHANES")]

GenusNHANES <- merge(GenusNHANES, x)
PhylumNHANES <- merge(PhylumNHANES, x)

rm("HealthyMicrobiome", "MyGenus",           "MyPhylum",          "nhanes",            "Pairs",            
   "PerioMicrobiome",        "psm_healthy",       "psm_perio",         "x" )

dput(names(PhylumNHANES))

names(GenusNHANES) <- c("SEQN", "Porphyromonas", "Prevotella", "Tannerella", "Aggregatibacter", 
                          "Fusobacterium", "Streptococcus", "Parvimonas", "Campylobacter", 
                          "Eikenella", "[Eubacterium] nodatum group", "Capnocytophaga", 
                          "Veillonella", "Actinomyces", "Prevotella_7", "Selenomonas", 
                          "Treponema", "Tax", "Run", "Sex", "AgeYrs", "Periodontitis")

names(PhylumNHANES) <-c("SEQN", "Bacteroidetes", "Proteobacteria", "Fusobacteriota", 
                        "Firmicutes", "Campylobacterota", "Bacteroidota", "Actinobacteriota", 
                        "Spirochaetota", "Tax", "Run", "Sex", "AgeYrs", 
                        "Periodontitis")

save.image("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/NHANES taxa.RData")


