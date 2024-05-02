# ================================================================================
#
# Retreaving SRAs for each BioProject included
#
# ================================================================================
library(rentrez)
library(XML)
library(dplyr)
library(readxl)
library(here)

# Setting NCBI API key
set_entrez_key("5fad35c7d12ae7d345787000e086295a4f08")

# Here the list of included BioProject IDs
my_project <- c("PRJNA552294", "PRJNA548383", "PRJNA508385")

# Store my metadata here:
MySRA <- data.frame()

# Performin search and retrieve metadata
retrieve_metadata <- function(project_id) {
  sra <- entrez_search(db = "sra", project_id, use_history = TRUE)
  cat(paste0("Retrieving ", sra$count, " SRAs information from BioProject ", project_id, "\n"))
  sra_xml <- entrez_fetch(db = "sra", web_history = sra$web_history, rettype = 'runinfo', retmode = 'xml', parsed = TRUE)
  sra_metadata <- xmlToDataFrame(sra_xml)
  return(sra_metadata)
}

# Applying the function to each project
MySRA <- bind_rows(lapply(my_project, retrieve_metadata))





# ================================================================================
#
# Setup directories and SRAs List to download
#
# ================================================================================
library(parallel)

# Setup directories and SRAs List to download
my_fastq <- MySRA$Run                                             # SRA list to download
sra_dir <- "/Users/davidepietropaoli/Downloads/WGS"               # Where SRA will be downloaded
fastq_out <- "/Users/davidepietropaoli/Downloads/FASTQ"           # Where FASTQ will be stored after SRA to FASTQ conversion
bowtie2out <- "/Users/davidepietropaoli/Downloads/bowtie2out"     # Where bowtie2out will be placed after Metaphlan4
MPA4 <- "/Users/davidepietropaoli/Downloads/MPA4"                 # Where Metaphlan final files will be stored
MPA4_path <- "/Users/davidepietropaoli/miniforge3/bin"            # Where "metaphlan" in stored: run into the terminal "which metaphlan"
nproc <- detectCores()-1                                          # Numbers of available cores. Leave 1 core free



# Function to run system commands with error handling
run_system_command <- function(command, message_prefix = "") {
  cat(paste0("\033[1;32m", message_prefix, command, "\033[0m\n"))
  result <- system(command)
  if (result != 0) {
    stop(paste("Error executing command:", command))
  }
}

# Function to classify FASTQ files
classify_fastq_files <- function(files) {
  single_end_files <- character(0)
  paired_end_files_1 <- character(0)
  paired_end_files_2 <- character(0)
  
  for (file in files) {
    if (grepl("_1.fastq$", file)) {
      paired_end_files_1 <- c(paired_end_files_1, file)
    } else if (grepl("_2.fastq$", file)) {
      paired_end_files_2 <- c(paired_end_files_2, file)
    } else if (grepl("\\.fastq$", file)) {
      single_end_files <- c(single_end_files, file)
    }
  }
  
  return(list(single_end = single_end_files, paired_end_1 = paired_end_files_1, paired_end_2 = paired_end_files_2))
}


# Caffeinate the system before running the prefetch command
system("caffeinate -i", wait = FALSE)  # -i: Prevent the system from idle sleeping.


# ================================================================================
# Download WGS sequences (SRA file) using the command "prefetch"
# ================================================================================
for (i in seq_along(my_fastq)) {
  run_system_command(paste("prefetch", my_fastq[i], "-p -v --verify yes --output-directory", file.path(sra_dir)), paste0("Downloading ", i, " of ", length(my_fastq), ": "))
  gc()
  Sys.sleep(1)
}

# ================================================================================
# Check SRA file integrity
# ================================================================================
for (i in seq_along(my_fastq)) {
  run_system_command(paste("vdb-validate", file.path(sra_dir, paste0(my_fastq[i], "/", my_fastq[i], ".sra"))), paste0("Checking integrity ", i, " of ", length(my_fastq), ": "))
  gc()
  }

# ================================================================================
# Extract SRA to FASTQ files
# ================================================================================
for (i in seq_along(my_fastq)) {
  run_system_command(paste("fasterq-dump -p", file.path(sra_dir, paste0(my_fastq[i], "/", my_fastq[i], ".sra")), "--outdir", fastq_out), paste0("Extracting FASTQ ", i, " of ", length(my_fastq), ": "))
}

# ================================================================================
# Verify if FASTQ files are paired-end or single-end
# ================================================================================
all_fastq <- list.files(fastq_out, pattern = "\\.fastq$", full.names = TRUE)
file_classification <- classify_fastq_files(all_fastq)

cat(paste0("Single-end FASTQ files: ",length(file_classification$single_end), "\n",
           "Paired-end R1 FASTQ files: ",length(file_classification$paired_end_1), "\n",
           "Paired-end R2 FASTQ files: ",length(file_classification$paired_end_2), "\n"))


# ================================================================================
# Performing MetaPhlAn4 for single-paired and for paired-end files
# ================================================================================
Sys.setenv(PATH = paste(MPA4_path, Sys.getenv("PATH"), sep = ":"))
metaphlan <- ("/Users/davidepietropaoli/miniforge3/bin/metaphlan ")

# Caffeinate the system before running the prefetch command
system("caffeinate -i", wait = FALSE)  # -i: Prevent the system from idle sleeping.

# Extract the desired part using regular expressions
matches <- regmatches(all_fastq, regexpr("SRR(\\d+)", all_fastq))
my_fastq <- unique(matches)

# Running MetaPhlAn4 for single end files
for (i in seq_along(file_classification$single_end)) {
  command <- paste0(metaphlan, file_classification$single_end[i],  
                    " --input_type fastq  --bowtie2out ",  bowtie2out, "/", my_fastq[i], ".bowtie2.bz2 -o ", 
                    MPA4, "/", my_fastq[i], "_metagenome.txt --nproc ", nproc)
  run_system_command(command)
}


# Running MetaPhlAn4 for pairend end files
for (i in seq_along(file_classification$paired_end_1)) {
  cat(paste0("Running MetaPhlAn: ", i, " of ", length(file_classification$paired_end_1), "\n"))
  command <- paste0(metaphlan, file_classification$paired_end_1[i], ",",file_classification$paired_end_2[i], 
                    " --input_type fastq  --bowtie2out ",  bowtie2out, "/", my_fastq[i], ".bowtie2.bz2 -o ", 
                    MPA4, "/", my_fastq[i], "_metagenome.txt --nproc ", nproc)
  run_system_command(command)
}


# ================================================================================
# Merging MetaPhlAn Tables 
# ================================================================================
command <- paste0("merge_metaphlan_tables.py ", MPA4, "/*_metagenome.txt > ", MPA4,"/MERGED_abundance_table.txt")
run_system_command(command)


# ================================================================================
# Convert MetaPhlAn Merged table to PhyloSeq object 
# ================================================================================
library(phyloseq)
library(devtools)
library(readr)
WGS_MetaData <- read_csv(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Retriving from PubMed/V2 WGS - MetaData.csv"), 
                            col_types = cols(...1 = col_skip()))

# Load my Library from GitHub
source_url("https://raw.githubusercontent.com/PietropaoliLab/r_functions/main/metaphlanToPhyloseq")


# Load MetaPhlan merged data and prepare it for conversion into PhyloSeq Obj
metaphlan <- read.table("~/Downloads/MPA4/MERGED_abundance_table.txt", header = TRUE, skip = 1)
rownames(metaphlan) <- metaphlan$clade_name     # Assign row names as clades
colnames(metaphlan) <- gsub("_metagenome", "", colnames(metaphlan))   # Standardize Colouns names and metadata names
metaphlan <- metaphlan[-1]  # Remove the first coloum

# Generating foo metadata dataset to prove the funcion
metadatadf <- as.data.frame(WGS_MetaData)
row.names(metadatadf) <- metadatadf$Run
sample <- metadatadf

wgs <- metaphlanToPhyloseq(metaphlan, metadat = sample)
wgs

plot_bar(wgs, fill = "Phylum")

save(wgs, file = here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS_for_revision.RData"))

load(here("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/WGS_for_revision.RData"))
