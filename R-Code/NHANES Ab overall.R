library(ggpubr)

# Load NHANES III
load("~/Google Drive/01- UNIVERSITA/AQ Projects/NHANES/NHANES III/data/NHANES for analysis DEC2018.RData")

myvars <- c("SEQN", "SEX", "AGE",  "periodontitis", "AT.LEAST.100.CIGS",
            "P.intermedia", "P.nigrescens", "T.forsythia", "A.actino.mix", "P.gingivalis",
            "F.nucleatum", "S.oralis",  "M.micros", "C.rectus", "E.corrodens", "E.nodatum", "S.intermedius", 
            "C.ochracea", "V.parvula", "A.naeslundii", "P.melaninogenica", 
            "S.noxia", "T.denticola", "S.mutans")

nhanes <- nhanes[myvars]

# Change names for consistency
names(nhanes)[2] <- "Sex"
names(nhanes)[3] <- "Age"
names(nhanes)[4] <- "Periodontitis"
names(nhanes)[5] <- "Smoking"


# Same factors for consistency
unique(nhanes$Sex)
nhanes$Sex <- as.character(nhanes$Sex)
nhanes$Sex <- gsub("Male", "male", nhanes$Sex)
nhanes$Sex <- gsub("Female", "female", nhanes$Sex)
nhanes$Periodontitis <- ifelse(nhanes$Periodontitis == 1, "Periodontitis", "Healthy")
nhanes$Smoking <- ifelse(nhanes$Smoking == 1, "Smokers", "Non-smokers")

# Excluding individuals without periodontitis status
nhanes <- subset(nhanes, !is.na(Periodontitis))

# Excluding individuals without smoking status
nhanes <- subset(nhanes, !is.na(Smoking))

# Excluding individuals without Age
nhanes <- subset(nhanes, !is.na(Age))

# Excluding individuals without Sex
nhanes <- subset(nhanes, !is.na(Sex))



# Performing Log10 normalization
MyStrains <- names(nhanes[6:24])
for (i in MyStrains) {
  foo <- log10(1 + nhanes[i]) # add +1 to avoid -Inf results
  names(foo) <- paste0(names(foo),".log10")
  nhanes <- cbind(nhanes, foo)
  rm(foo)
}
rm(MyStrains, i)







#######################################
#
## Generating Genus Ab in NHANES
#
#######################################

# Eubacterium = E.nodatum
nhanes$`[Eubacterium] nodatum group` <- nhanes$E.nodatum.log10

# Actinomyces = A.naeslundii
nhanes$Actinomyces <- nhanes$A.naeslundii.log10

# Aggregatibacter = a.actimo.mix
nhanes$Aggregatibacter <- nhanes$A.actino.mix.log10

# Campylobacter = C.rectus
nhanes$Campylobacter <- nhanes$C.rectus.log10

# Capnocytophaga = C.ocracea
nhanes$Capnocytophaga <- nhanes$C.ochracea.log10

# Eikenella = E.corrodens
nhanes$Eikenella <- nhanes$E.corrodens.log10

# Fusobacterium = F.nucleatum
nhanes$Fusobacterium <- nhanes$F.nucleatum.log10

# Parvimonas = M.micros
nhanes$Parvimonas <- nhanes$M.micros.log10

# Porphyromonas = P.gingivalis
nhanes$Porphyromonas <- nhanes$P.gingivalis.log10

# Prevotella = P.intermedia; P.nigrescens;
nhanes$Prevotella <- rowSums(nhanes[,c("P.intermedia.log10", "P.nigrescens.log10")])

# Prevotella_7 = P.melaninogenica
nhanes$Prevotella_7 <- nhanes$P.melaninogenica.log10

# Selenomonas = S.noxia
nhanes$Selenomonas <- nhanes$S.noxia.log10

# Streptococcus = S.oralis; S.intermedius; S.mutans
nhanes$Streptococcus <- rowSums(nhanes[,c("S.oralis.log10", "S.intermedius.log10", "S.mutans.log10")])

# Tannerella = T.forsythia
nhanes$Tannerella <- nhanes$T.forsythia.log10

# Treponema = T.denticola
nhanes$Treponema <- nhanes$T.denticola.log10

# Veillonella = V.parvula
nhanes$Veillonella <- nhanes$V.parvula.log10


#######################################
#
## Generating Phylum Ab in NHANES
#

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




# Subsetting TAXA according to Phylum in NHANES antibodies
MyPhylum <- c("Bacteroidetes", "Proteobacteria", "Fusobacteriota", "Firmicutes", 
              "Campylobacterota", "Bacteroidota", "Actinobacteriota", "Spirochaetota")

MyGenus <- c("Porphyromonas", "Prevotella", "Tannerella", "Aggregatibacter", 
             "Fusobacterium", "Streptococcus", "Parvimonas", "Campylobacter", 
             "Eikenella", "[Eubacterium] nodatum group", "Capnocytophaga", 
             "Veillonella", "Actinomyces", "Prevotella_7", "Selenomonas", 
             "Treponema")





# Selecting only TAX of interest
PhylumNHANES <- nhanes[,c("SEQN","Sex", "Age", "Periodontitis", "Smoking", MyPhylum)]
PhylumNHANES$Tax <-"Phylum"

GenusNHANES <-  nhanes[,c("SEQN","Sex", "Age", "Periodontitis", "Smoking", MyGenus)]
GenusNHANES$Tax <-"Genus"


library(reshape2)
PhylumNHANES <- melt(PhylumNHANES, id.vars = c("SEQN", "Sex", "Age", "Periodontitis", "Smoking","Tax"))
GenusNHANES <- melt(GenusNHANES, id.vars = c("SEQN", "Sex", "Age", "Periodontitis", "Smoking","Tax"))




library(ggpubr)
MyTheme <-   theme_bw(10)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 2/1,
        legend.key = element_blank(),
        legend.background=element_blank())



ggboxplot(GenusNHANES, x="variable", y="value",
          color = "Sex", palette = c("#F88B9D", "#8ECEFD"),
          add = "jitter",
          facet.by = c("Periodontitis", "Smoking"),
          add.params = list(size = 1, alpha = .3))+
  stat_compare_means(aes(color = Sex), method = "wilcox.test", label = "p.signif", show.legend = FALSE, label.y = 10.4)+
  labs(x="", y="Ab titers (Log10)", title = "Ab - Genus")+
  #scale_y_continuous(trans = "log10")+
  coord_cartesian(ylim = c(0.1, 10.6))+
  MyTheme+
  theme(aspect.ratio = 1/3,
        axis.text.x = element_text(size=8, angle = 35, hjust = 1))


GenusNHANES$AgeCat <- cut(GenusNHANES$Age, breaks = c(0, 45, 50, 55, 60, 65, 70, Inf), labels = c("<45", "<50", "<55", "<60", "<65", "<70", "70+"))
  
ggline(GenusNHANES, x = "AgeCat", y = "value", add = "mean_se", color = "Sex",
       facet.by = c("Periodontitis", "variable"))+
  MyTheme

