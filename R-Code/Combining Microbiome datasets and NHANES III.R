library(phyloseq)
# Load my 16S data
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Load NHANES III
load("AQ Projects/NHANES/NHANES III/data/NHANES for analysis DEC2018.RData")

# Working on microbiome data
# PhyloSeq metadata
MetaDataHealthy <- data.frame(psm_healthy@sam_data)
MetaDataPerio <- data.frame(psm_perio@sam_data)
MetaDataHealthy <- MetaDataHealthy[, c("Run", "Sex", "AgeYrs", "Periodontitis", "Site", "Smoking")]
MetaDataPerio <- MetaDataPerio[, c("Run", "Sex", "AgeYrs", "Periodontitis", "Site", "Smoking")]
MetaData <- rbind(MetaDataHealthy, MetaDataPerio)

# Select only individuals with Subgingival microbiome
# Before subsetting 570 individuals
# After subsetting 199 
MetaData <- subset(MetaData, Site == "Subgingival")

# After subsetting unknown smoking habits
# we have 105
MetaData <- subset(MetaData, !is.na(Smoking))

MetaData <- MetaData[-5]

# Screening data
table(MetaData$Periodontitis, MetaData$Sex)

# some cleaning
rm(MetaDataHealthy, MetaDataPerio)


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
Healthy <- subset(Pairs, Periodontitis_16S == "Healthy")

# Periodontitis pairs
MyRun <- Perio$Run
PerioMicrobiome <- subset_samples(psm_perio, Run %in% MyRun)


# Healthy pairs
MyRun <- Healthy$Run
HealthyMicrobiome <- subset_samples(psm_healthy, Run %in% MyRun)


# NHANES pairs
# Load NHANES III
load("AQ Projects/NHANES/NHANES III/data/NHANES for analysis.RData")
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

x <- Pairs[,c("Run", "SEQN")]

GenusNHANES <- merge(GenusNHANES, x)
PhylumNHANES <- merge(PhylumNHANES, x)



library(reshape2)
PhylumNHANES <- melt(PhylumNHANES, id.vars = c("SEQN", "Run", "Tax"))
GenusNHANES <- melt(GenusNHANES, id.vars = c("SEQN", "Run", "Tax"))

names(PhylumNHANES) <- c("SEQN","Run","Tax","Phylum", "value" )  
names(GenusNHANES) <- c("SEQN","Run","Tax","Genus", "value" )  

PhylumNHANES$Confirm <- PhylumNHANES$Phylum
GenusNHANES$Confirm <- GenusNHANES$Genus


# Working on 16S datasets
HealthyMicrobiomeGenus <- subset_taxa(HealthyMicrobiome, Genus %in% MyGenus)
PerioMicrobiomeGenus <- subset_taxa(PerioMicrobiome, Genus %in% MyGenus)

HealthyMicrobiomePhylum <- subset_taxa(HealthyMicrobiome, Phylum %in% MyPhylum)
PerioMicrobiomePhylum <- subset_taxa(PerioMicrobiome, Phylum%in% MyPhylum)


# Convert in relative abundance. "Compositional"
library(microbiome)
HealthyMicrobiomeGenus <- transform(HealthyMicrobiomeGenus, transform = "compositional")
PerioMicrobiomeGenus <- transform(PerioMicrobiomeGenus, transform = "compositional")

HealthyMicrobiomePhylum <- transform(HealthyMicrobiomePhylum, transform = "compositional")
PerioMicrobiomePhylum <- transform(PerioMicrobiomePhylum, transform = "compositional")

HealthyGenus <- psmelt(HealthyMicrobiomeGenus)
PerioGenus <- psmelt(PerioMicrobiomeGenus)

HealthyPhylum <- psmelt(HealthyMicrobiomePhylum)
PerioPhylum <- psmelt(PerioMicrobiomePhylum)

MyVars <- c("Run", "Sex", "AgeYrs", "Periodontitis", "Site", "Smoking", "Abundance",
            "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

HealthyGenus <- HealthyGenus[MyVars]
PerioGenus <- PerioGenus[MyVars]
HealthyPhylum <- HealthyPhylum[MyVars]
PerioPhylum <- PerioPhylum[MyVars]


CombinedGenusPeriodontitis <- left_join(PerioGenus, GenusNHANES, by =c("Run", "Genus"))
CombinedGenusHealthy <- left_join(HealthyGenus, GenusNHANES, by =c("Run", "Genus"))


CombinedPhylumPeriodontitis <- left_join(PerioPhylum, PhylumNHANES, by =c("Run", "Phylum"))
CombinedPhylumHealthy <- left_join(HealthyPhylum, PhylumNHANES, by =c("Run", "Phylum"))

names(CombinedGenusPeriodontitis)[names(CombinedGenusPeriodontitis) == 'value'] <- "Ab_titers"
names(CombinedGenusHealthy)[names(CombinedGenusHealthy) == 'value'] <- "Ab_titers"
names(CombinedPhylumPeriodontitis)[names(CombinedPhylumPeriodontitis) == 'value'] <- "Ab_titers"
names(CombinedPhylumHealthy)[names(CombinedPhylumHealthy) == 'value'] <- "Ab_titers"

# Define smokers and non smokers
CombinedGenusHealthy$Smoking <- ifelse(CombinedGenusHealthy$Smoking == "No", "Non-smokers", "Smokers")
CombinedGenusPeriodontitis$Smoking <- ifelse(CombinedGenusPeriodontitis$Smoking == "No", "Non-smokers", "Smokers")

library(ggpubr)

MyTheme <-   theme_bw(10)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 1/1,
        legend.key = element_blank(),
        legend.background=element_blank())



p1 <- ggplot(CombinedGenusHealthy, aes(x=Ab_titers, y=Abundance, group = Sex, color = Sex))+
  geom_point(color = "#d9d9d9", size = 2)+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  geom_smooth(aes(color = Sex), method = "lm", se = FALSE)+
  # facet_wrap(~ Smoking)+
  xlim(0, 0.45)+ 
  ylim(0, 1)+
  MyTheme+
  labs(x="Ab titers (%)", y="Genus relative abundance (%)", title = "Healhy", subtitle = "Genus")+
  stat_cor(aes(color = Sex), method = "pearson", label.x.npc =.05, label.y.npc  = .95, show.legend = FALSE, size = 3)+
  annotate("text", x = 0.35, y = .15, size = 3, label = paste0("n=", length(unique(CombinedGenusHealthy$SEQN))))


p2 <- ggplot(CombinedGenusPeriodontitis, aes(x=Ab_titers, y=Abundance, group = Sex, color = Sex))+
  geom_point(color = "#d9d9d9", size = 2)+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  geom_smooth(aes(color = Sex), method = "lm", se = FALSE)+
  # facet_wrap(~ Smoking)+
  xlim(0, 0.45)+ 
  ylim(0, 1)+
  MyTheme+
  labs(x="Ab titers (%)", y="Genus relative abundance (%)", title = "Periodontitis", subtitle = "Genus")+
  stat_cor(aes(color = Sex), method = "pearson", label.x.npc =.05, label.y.npc  = .95, show.legend = FALSE, size = 3)+
  annotate("text", x = 0.35, y = .15, size = 3, label = paste0("n=", length(unique(CombinedGenusPeriodontitis$SEQN))))



p3 <- ggplot(CombinedPhylumHealthy, aes(x=Ab_titers, y=Abundance, group = Sex, color = Sex))+
  geom_point(color = "#d9d9d9", size = 2)+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  geom_smooth(aes(color = Sex), method = "lm", se = FALSE)+
  # facet_wrap(~ Smoking)+
  xlim(0, 0.45)+ 
  ylim(0, 1)+
  MyTheme+
  labs(x="Ab titers (%)", y="Phylum relative abundance (%)", title = "Healhy", subtitle = "Phylum")+
  stat_cor(aes(color = Sex), method = "pearson", label.x.npc =.05, label.y.npc  = .95, show.legend = FALSE, size = 3)+
  annotate("text", x = 0.35, y = .15, size = 3, label = paste0("n=", length(unique(CombinedPhylumHealthy$SEQN))))


p4 <- ggplot(CombinedPhylumPeriodontitis, aes(x=Ab_titers, y=Abundance, group = Sex, color = Sex))+
  geom_point(color = "#d9d9d9", size = 2)+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  geom_smooth(aes(color = Sex), method = "lm", se = FALSE)+
  # facet_wrap(~ Smoking)+
  xlim(0, 0.45)+ 
  ylim(0, 1)+
  MyTheme+
  labs(x="Ab titers (%)", y="Phylum relative abundance (%)", title = "Periodontitis", subtitle = "Phylum")+
  stat_cor(aes(color = Sex), method = "pearson", label.x.npc =.05, label.y.npc  = .95, show.legend = FALSE, size = 3)+
  annotate("text", x = 0.35, y = .15, size = 3, label = paste0("n=", length(unique(CombinedPhylumPeriodontitis$SEQN))))


# pdf(file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Scatter Plot Relative Abundance vs Ab titer NHANES.pdf", width = 7, height = 7)


library(patchwork)
p1+p2+p3+p4+ plot_layout(guides = "collect")

  # ggexport(filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Scatter Plot Relative Abundance vs Ab titer NHANES.pdf")






ggboxplot(CombinedGenusHealthy, x="Genus", y="Ab_titers",
          color = "Sex", palette = c("#F88B9D", "#8ECEFD"),
          add = "jitter",
          add.params = list(size = 1.5, alpha = .3))+
  stat_compare_means(aes(color = Sex), method = "wilcox.test", label = "p.signif", show.legend = FALSE)+
  labs(x="", y="Ab titers (%)", title = "Healthy", subtitle = "Genus")+
  MyTheme+
  theme(aspect.ratio = 1/3,
        axis.text.x = element_text(size=8, angle = 35, hjust = 1))



ggboxplot(CombinedGenusPeriodontitis, x="Genus", y="Ab_titers",
          color = "Sex", palette = c("#F88B9D", "#8ECEFD"),
          add = "jitter",
          add.params = list(size = 1.5, alpha = .3))+
  stat_compare_means(aes(color = Sex), method = "wilcox.test", label = "p.signif", show.legend = FALSE)+
  labs(x="", y="Ab titers (%)", title = "Periodontitis", subtitle = "Genus")+
  MyTheme+
  theme(aspect.ratio = 1/3,
        axis.text.x = element_text(size=8, angle = 35, hjust = 1))



ggboxplot(CombinedPhylumHealthy, x="Phylum", y="Ab_titers",
          color = "Sex", palette = c("#F88B9D", "#8ECEFD"), 
          add = "jitter",
          add.params = list(size = 1.5, alpha = .3))+
  stat_compare_means(aes(color = Sex), method = "wilcox.test", label = "p.signif", show.legend = FALSE)+
  labs(x="", y="Ab titers (%)", title = "Healthy", subtitle = "Phylum")+
  MyTheme+
  theme(aspect.ratio = 1/3,
        axis.text.x = element_text(size=8, angle = 35, hjust = 1))



ggboxplot(CombinedPhylumPeriodontitis, x="Phylum", y="Ab_titers",
          color = "Sex", palette = c("#F88B9D", "#8ECEFD"), 
          add = "jitter",
          add.params = list(size = 1.5, alpha = .3))+
  stat_compare_means(aes(color = Sex), method = "wilcox.test", label = "p.signif", show.legend = FALSE)+
  labs(x="", y="Ab titers (%)", title = "Periodontitis", subtitle = "Phylum")+
  MyTheme+
  theme(aspect.ratio = 1/3,
        axis.text.x = element_text(size=8, angle = 35, hjust = 1))




