library(phyloseq)
library(ggpubr)
library(ggplot2)
library(microbiome)
library(plyr)
library(reshape2)
library(microbiomeMarker)
library(extrafont)
library(showtext)
library(rentrez)
library(XML)
library(dplyr)
library(readxl)
library(here)
loadfonts(device = "pdf")

# Load my 16S data
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")


my_project <- c("PRJNA774299", "PRJNA774981")
p <- subset_samples(psm_perio, BioProject %in% my_project)

######################
#
# Retrieving metadata for each subject
#
######################

# Remove duplicates based on BioSamples...N=10301 # BioSamples duplicated = 317
MyBioSamples <- p@sam_data$BioSample


# Counter
n <- length(MyBioSamples)

# Defining dataframe
MetaData <- data.frame()

# Inizialize the cycle counter
k = 1

# Take system time 
stime <- Sys.time()

for (i in MyBioSamples) {
  
  itime <- Sys.time()
  cat(paste0("Cycle ", k, " - Retrieving metadata from Subject: ", i, "\n"))  
  
  # Incrementing counter
  k = k + 1
  # i = "SAMN12211471"
  sra <- entrez_search(db="biosample", i, use_history=TRUE)
  
  # Storing my data info XML file
  sra_xml <- entrez_fetch(db = "biosample", web_history = sra$web_history, rettype='xml', parsed=TRUE)
  
  
  # Extracting name of variables metadata
  AvailableMetaData <-  unlist(lapply(xpathApply(sra_xml, '//Attribute', xmlAttrs), '[[', "attribute_name"))
  
  # Extracting metadata for each Subject
  SubjectMetaData <- xpathSApply(sra_xml, "//Attribute[@attribute_name]", xmlValue)
  
  foo <- paste(AvailableMetaData, SubjectMetaData, sep = ": ")
  foo <- paste(foo, collapse = " | ")
  
  tmp <- as.data.frame(cbind("BioSample" = i))
  tmp <- cbind(tmp, foo)
  
  MetaData <- rbind(MetaData, tmp)
  
  
  # Add some  dalay in order to retreive max 3 tables in a second - Otherwise NCBI gives error
  # Sys.sleep(0.06) 
  
  # get partial time
  ptime <- Sys.time()
  
  # Print time needed to retrieve project info
  cat(paste0("Data retrieved in: ", round(difftime(ptime, itime, unit = "sec"), 2), " sec.  -  ", round(((k/n)*100), 1),"% Done.\n\n"))
  cat(paste0("\n"))
  Sys.sleep(0.1)
}

# Print time needed to retrieve ALL projects info
etime <- Sys.time()
cat(paste0("OPERATION CONCLUDED IN: ", round(difftime(etime, stime, unit = "min"), 2), " minutes","\n"))

rm(stime, ptime, etime, itime)

# Defining names
names(MetaData) <- c("BioSample", "MetaData")
rownames(MetaData) <- MetaData$BioSample

#--------------------------
## Starting producing dataframe of metadata
# Appending BioSamples ID
Biosample <- paste0("BioSample: ", MetaData$BioSample, " | ")

# Preparing for pasting Biosample with MetaData
TextMetaData <- MetaData$MetaData

# Pasting Biosample with MetaData
TextMetaDataWithBioSample <- paste0(Biosample, TextMetaData)

# Cleaning
rm(Biosample, TextMetaData)

# Converting into dataset
SubjectsMetaData <- type.convert(as.data.frame(read.dcf(
  textConnection(paste(gsub("\\s+\\|\\s+", "\n", TextMetaDataWithBioSample), 
                       collapse="\n\n")))), as.is = TRUE)


MetaData <- SubjectsMetaData[,c(1, 14, 15, 17)]
old_metadata <- data.frame(sample_data(p))

# Now you can attempt to merge the sample data again
merged_data <- merge(old_metadata, MetaData, by="BioSample", all.x=TRUE)
rownames(merged_data) <- merged_data$Run
sample_data(p) <- merged_data

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
load("~/Il mio Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Periodontitis_according_to_severity.RData")

# Defining my theme
MyTheme <-  theme_bw(base_size = 11)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 2.5/1,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir"))

table(p@sam_data$Periodontitis_extent, p@sam_data$Sex)
p@sam_data$Sex <- factor(p@sam_data$Sex, levels = c("female", "male"), labels = c("F", "M"))

# Plot periodontitis
p1 <- plot_richness(p, x="Sex", measures="Observed", color = "Sex")+
  facet_grid(~ Periodontitis_extent, scales = "free")+  # Adjust facet_wrap to use scales = "free"
  geom_violin(aes(fill=Sex), alpha=.3)+
  geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=1.3, stroke =.4)+
  geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
  labs(title = "A", y= "Observed alpha diversity")+
  #coord_fixed(ylim = c(0,4))+
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme+
  stat_compare_means(aes(group = Periodontitis_extent), comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family = "Avenir")

p2 <- plot_richness(p, x="Sex", measures="Shannon", color = "Sex")+
  facet_grid(~ Periodontitis_extent, scales = "free")+  # Adjust facet_wrap to use scales = "free"
  geom_violin(aes(fill=Sex), alpha=.3)+
  geom_point(position = position_jitter(seed = 101, width = .20), alpha = I(.7), shape=21, fill="gray90", size=1.3, stroke =.4)+
  geom_boxplot(width =.15, fill = "white", outlier.shape = NA)+
  labs(title = "B", y= "Shannon diversity")+
  #coord_fixed(ylim = c(0,4))+
  scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
  scale_color_manual(values = c("#F88B9D", "#8ECEFD"))+
  MyTheme+
  stat_compare_means(aes(group = Periodontitis_extent), comparisons = list(c("F", "M")), method = "wilcox.test", label = "p.signif", family = "Avenir")


p1 + p2







