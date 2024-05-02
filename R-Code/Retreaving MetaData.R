# Load libraries
library(rentrez)
library(XML)
library(dplyr)
library(purrr)

# Setting NCBI API key
set_entrez_key("5fad35c7d12ae7d345787000e086295a4f08")

# Function to retrieve metadata for a given project ID
retrieve_metadata <- function(project_id) {
  cat(paste0("Retrieving subject IDs from project: ", project_id, "\n"))
  
  # Performing search
  sra <- entrez_search(db = "sra", project_id, use_history = TRUE)
  
  # Storing data info in XML file
  sra_xml <- entrez_fetch(db = "sra", web_history = sra$web_history, rettype = 'runinfo', retmode = 'xml', parsed = TRUE)
  
  # Convert XML to DataFrame
  sra_metadata <- xmlToDataFrame(sra_xml)[, c("Run", "BioProject", "BioSample")]
  
  return(sra_metadata)
}

# Retrieve metadata using purrr and dplyr
BioSamples <- map_df(my_sra, retrieve_metadata)

# Print the final result
print(BioSamples)

# Define function to retrieve metadata from a BioSample ID
retrieve_biosample_metadata <- function(biosample_id) {
  cat(paste0("Retrieving metadata from BioSample: ", biosample_id, "\n"))
  
  sra <- entrez_search(db = "biosample", biosample_id, use_history = TRUE)
  
  # Storing data info XML file
  sra_xml <- entrez_fetch(db = "biosample", web_history = sra$web_history, rettype='xml', parsed=TRUE)
  
  # Extracting name of variables metadata
  available_metadata <-  unlist(lapply(xpathApply(sra_xml, '//Attribute', xmlAttrs), '[[', "attribute_name"))
  
  # Extracting metadata for each BioSample
  subject_metadata <- xpathSApply(sra_xml, "//Attribute[@attribute_name]", xmlValue)
  
  metadata_text <- paste(available_metadata, subject_metadata, sep = ": ")
  metadata_text <- paste(metadata_text, collapse = " | ")
  
  return(metadata_text)
}

# Retrieve metadata for each BioSample
metadata_list <- map_chr(BioSamples$BioSample, retrieve_biosample_metadata)

# Combine metadata with BioSamples
SubjectsMetaData <- cbind(BioSamples, MetaData = type.convert(as.data.frame(read.dcf(
                                                            textConnection(paste(gsub("\\s+\\|\\s+", "\n", metadata_list), 
                                                            collapse="\n\n")))), as.is = TRUE))

# Print the final result
print(SubjectsMetaData)
