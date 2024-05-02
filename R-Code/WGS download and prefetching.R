# ================================================================================
#
# Setup directories and SRAs List to download
#
# ================================================================================
sra_dir <- "/Users/davidepietropaoli/Downloads/WGS/"              # My SRA input destination folder
fastq_out <- "/Users/davidepietropaoli/Downloads/FASTQ/"          # My FASTQ output destination folder
my_fastq <- c("SRR9641934", "SRR9641941")                         # My SRA to download
bowtie2out <- "/Users/davidepietropaoli/Downloads/bowtie2out/"    # bowtie2out files
MPA4 <- "/Users/davidepietropaoli/Downloads/MPA4/"                # Metaphlan 4 generated files



# ================================================================================
#
# Download my WGS sequences (SRA file) using the command "prefetch"
#
# ================================================================================
# Cycle for downloading my fastq files 
i = 1
for (fastq in my_fastq) {
  
  # Caffeinate the system before running the prefetch command
  system("caffeinate -i", wait = FALSE)  # -i: Prevent the system from idle sleeping.
  
  # Print the command
  command <- paste("prefetch", fastq, "-p -v --verify yes  --output-directory", 
                   paste0(sra_dir, fastq), sep = " ")
  cat("\033[1;32m", paste0(command, " - ", i, " of ",   length(my_fastq)), "\033[0m\n") # \033[1;32m sets the color to green; \033[0m resets the color to the default.
  
  # Run the prefetch command
  system(command)
  i = i+1
  
}

# ================================================================================
#
# Check SRA file integrity
#
# ================================================================================
# Check SRA file integrity
i = 1
for (fastq in my_fastq) {
  # Caffeinate the system before running the prefetch command
  system("caffeinate -i", wait = FALSE)  # -i: Prevent the system from idle sleeping.
  # Print the command
  command <- paste0("vdb-validate ",sra_dir, fastq, "/", fastq, ".sra") 
  cat("\033[1;32m", paste0(command, " - ", i, " of ",   length(my_fastq)), "\033[0m\n") # \033[1;32m sets the color to green; \033[0m resets the color to the default.
  
  # Run the vdb-validate command
  system(command)
  i = i+1
  
}


# ================================================================================
#
# Extracts SRA in FASTQ files 
#
# ================================================================================
# Cycle for EXTRACTING my fastq files 
i = 1
for (fastq in my_fastq) {
  
  # Caffeinate the system before running the prefetch command
  system("caffeinate -i", wait = FALSE)  # -i: Prevent the system from idle sleeping.

  # Print the command
  command <- paste0("fasterq-dump -p ",sra_dir, fastq, "/", fastq, ".sra", " --outdir ", fastq_out) 
  cat("\033[1;32m", paste0(command, " - ", i, " of ",   length(my_fastq)), "\033[0m\n") # \033[1;32m sets the color to green; \033[0m resets the color to the default.
  
  # EXTRACTING FASTQ
  system(command)
  i = i+1
  
}

# ================================================================================
#
#   Verify if FASTQ files are paired-end or single-end
#
# ================================================================================
fastq_out <- "/Users/davidepietropaoli/Downloads/FASTQ/"  # Corrected FASTQ output destination folder

# List all files in the folder with .fastq extension
all_fastq <- list.files(fastq_out, pattern = "\\.fastq$", full.names = TRUE)

# Initialize vectors for single and double paired fastq files
single_end_files <- character(0)
paired_end_files_1 <- character(0)
paired_end_files_2 <- character(0)

# Loop through all files and classify them
for (file in all_fastq) {
  file <- file.path(fastq_out, file)  # Use file.path to concatenate paths
  if (grepl("_1.fastq$", file)) {
    paired_end_files_1 <- c(paired_end_files_1, file)
  } else if (grepl("_2.fastq$", file)) {
    paired_end_files_2 <- c(paired_end_files_2, file)
  } else if (grepl("\\.fastq$", file)) {
    single_end_files <- c(single_end_files, file)
  }
}

# Print or use the vectors as needed
print("Single-end files:")
print(single_end_files)

print("Paired-end files (R1):")
print(paired_end_files_1)

print("Paired-end files (R2):")
print(paired_end_files_2)




# TO ACTIVATE THE ENVIRONMENT FOR METAPHLAN4
# " cd ~"
# "source mpa4/bin/activate"



fastq <- "SRR9641934"
Sys.setenv(PATH = paste("/Users/davidepietropaoli/miniforge3/bin", Sys.getenv("PATH"), sep = ":"))
metaphlan <- ("/Users/davidepietropaoli/miniforge3/bin/metaphlan ")
# metaphlan /Users/davidepietropaoli/Downloads/FASTQ/SRR9641934_1.fastq,/Users/davidepietropaoli/Downloads/FASTQ/SRR9641934_2.fastq  --input_type fastq  --bowtie2out /Users/davidepietropaoli/Downloads/bowtie2out/metagenome.bowtie2.bz2 -o /Users/davidepietropaoli/Downloads/MPA4/profiled_metagenome.txt --nproc 8
command <- paste0(metaphlan, fastq_out, "/",fastq, "_1.fastq,",  fastq_out, "/", fastq, "_2.fastq ",  "--input_type fastq  --bowtie2out ",  bowtie2out, "/", fastq, ".bowtie2.bz2 -o ", MPA4, "/", fastq, "_metagenome.txt --nproc 8")



command
system(command)

system("/Users/davidepietropaoli/miniforge3/bin/metaphlan --version")






system(command)
