tax <- data.frame(tax_table(ps))

library(stringr)
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean$kiss <- NA

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fille holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  
  ## Fill in missing taxonomy
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax.clean <- tax.clean[,c(1:6,8)]               


n_seqs <- seq(1, nrow(tax.clean))
len_n_seqs <- nchar(max(n_seqs))

taxa_names <- paste("OTU", formatC(n_seqs, 
                                       width = len_n_seqs, 
                                       flag = "0"), sep = "_")

taxa_names <- paste0(taxa_names, "_", tax.clean$Species)

rownames(tax.clean) <- taxa_names
tax.clean <- tax.clean[1:6]

tax.clean <- as.matrix(tax.clean)

taxa_names(ps) <- taxa_names

head(taxa_names(ps))



