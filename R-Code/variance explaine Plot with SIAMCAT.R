library(tidyverse)
library(SIAMCAT)
library(ggplot2)
library(phyloseq)
library(microbiome)

load('G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData')
ps <- psm_perio
sample_data(ps) <- sample_data(ps)[,c('BioProject', 'Sex','Site', 'Smoking', 'size_MB')]
sample_data(ps)$size_MB <- as.numeric(sample_data(ps)$size_MB)



ps <- phyloseq::tax_glom(ps, taxrank = "Genus")
ps <- microbiome::transform(ps, "compositional")

label <- create.label(meta=sample_data(ps),
                      label = "Sex",
                      case = "female")

# run the constructor function
sc.obj <- siamcat(phyloseq=ps, label=label)

check.confounders(sc.obj, fn.plot = 'G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Confounder_plot_Perio_meta.pdf',
                  feature.type='original')


###################################################################################



load('G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData')
ps <- psm_healthy
sample_data(ps) <- sample_data(ps)[,c('BioProject', 'Sex','Site', 'Smoking', 'size_MB')]
sample_data(ps)$size_MB <- as.numeric(sample_data(ps)$size_MB)



ps <- phyloseq::tax_glom(ps, taxrank = "Genus")
ps <- microbiome::transform(ps, "compositional")

label <- create.label(meta=sample_data(ps),
                      label = "Sex",
                      case = "female")

# run the constructor function
sc.obj <- siamcat(phyloseq=ps, label=label)

check.confounders(sc.obj, fn.plot = 'G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Confounder_plot_Healthy_meta.pdf',
                  feature.type='original')











