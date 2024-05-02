# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
library(rentrez)
library(XML)
library(readxl)
library(tidyverse)
library(here)


# Setting NCBI API key
set_entrez_key("5fad35c7d12ae7d345787000e086295a4f08")


######################
#
# Retrieving metadata for each subject
#
######################

# Remove duplicates based on BioSamples...N=10301 # BioSamples duplicated = 317
MyBioSamples <- MySRA$BioSample


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
  
  # Extacting Project ID
  ProjectID <- MySRA[which(MySRA$BioSample==i),]["BioProject"]
  
  foo <- paste(AvailableMetaData, SubjectMetaData, sep = ": ")
  foo <- paste(foo, collapse = " | ")
  
  tmp <- as.data.frame(cbind(ProjectID, "BioSample" = i))
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
names(MetaData) <- c("BioProject", "BioSample", "MetaData")
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


# Saving Metadata on CSV file
write.csv(SubjectsMetaData, here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Retriving from PubMed/V2 WGS - MetaData.csv"))


MyColumns <- c("BioSample", "isolation_source", "geo_loc_name",
               "ethnicity", "host_sex", "host_subject_id", "host_age")


SubjectsMetaData <- SubjectsMetaData[MyColumns]
SubjectsMetaData$Disease <- ifelse(grepl("^GAP|^LAP", SubjectsMetaData$host_subject_id), "Periodontitis", 
                              ifelse(grepl("^HS|^HNS", SubjectsMetaData$host_subject_id), "Healthy", "Other"))


SubjectsMetaData$Smoke <- ifelse(grepl("^HS", SubjectsMetaData$host_subject_id), "Smoker", 
                                   ifelse(grepl("^HNS", SubjectsMetaData$host_subject_id), "Non Smoker", NA))


SubjectsMetaData <- SubjectsMetaData[c(1:5, 7:9)]
dput(names(SubjectsMetaData))
names(SubjectsMetaData) <- c("BioSample", "Site", "Geo", "Eth", "Sex", "Age", "Disease", "Smoke")
SubjectsMetaData$Site <- ifelse(SubjectsMetaData$Site == "subgingival plaque", "Subgingival", SubjectsMetaData$Site )

MetaData <- merge(MySRA, SubjectsMetaData, by = "BioSample")

# Saving Metadata on CSV file
write.csv(MetaData, here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Retriving from PubMed/V2 WGS - MetaData.csv"))

