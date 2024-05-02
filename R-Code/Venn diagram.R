# Load library
library(ggVennDiagram)
library(xlsx)
library(reshape2)
library(tidyverse)
library(plyr)
library(extrafont)
library(showtext)
library(patchwork)
loadfonts(device = "pdf")

Sankey <- read.xlsx(file = "AQ Projects/Microbiome/Gender features of periodontal microbiome/Sankey data.xlsx", sheetIndex = 1)
Sankey <- Sankey[-1]

Sankey <- Sankey[,c(1,3,4,2)]

# Perio Smokers
if (length(subset(Sankey, Site == "Saliva" & Smoking == "Smokers" & Group == "Periodontitis")$Genus) <1 ) {
  Saliva = ""
}  else {
  Saliva = subset(Sankey, Site == "Saliva" & Smoking == "Smokers" & Group == "Periodontitis")$Genus
}

if (length(subset(Sankey, Site == "Plaque" & Smoking == "Smokers" & Group == "Periodontitis")$Genus) <1 ) {
  Plaque = ""
}  else {
  Plaque = subset(Sankey, Site == "Plaque" & Smoking == "Smokers" & Group == "Periodontitis")$Genus
}

if (length(subset(Sankey, Site == "Subgingival" & Smoking == "Smokers" & Group == "Periodontitis")$Genus) <1 ) {
  Subgingival = ""
}  else {
  Subgingival = subset(Sankey, Site == "Subgingival" & Smoking == "Smokers" & Group == "Periodontitis")$Genus
}
PS_Saliva <- Saliva
PS_Plaque <- Plaque
PS_Subgingival <- Subgingival
rm(Saliva, Plaque, Subgingival)



# Perio Non Smokers
if (length(subset(Sankey, Site == "Saliva" & Smoking == "Non Smokers" & Group == "Periodontitis")$Genus) <1 ) {
  Saliva = ""
}  else {
  Saliva = subset(Sankey, Site == "Saliva" & Smoking == "Non Smokers" & Group == "Periodontitis")$Genus
}

if (length(subset(Sankey, Site == "Plaque" & Smoking == "Non Smokers" & Group == "Periodontitis")$Genus) <1 ) {
  Plaque = ""
}  else {
  Plaque = subset(Sankey, Site == "Plaque" & Smoking == "Non Smokers" & Group == "Periodontitis")$Genus
}

if (length(subset(Sankey, Site == "Subgingival" & Smoking == "Non Smokers" & Group == "Periodontitis")$Genus) <1 ) {
  Subgingival = ""
}  else {
  Subgingival = subset(Sankey, Site == "Subgingival" & Smoking == "Non Smokers" & Group == "Periodontitis")$Genus
}

PNS_Saliva <- Saliva
PNS_Plaque <- Plaque
PNS_Subgingival <- Subgingival
rm(Saliva, Plaque, Subgingival)


# Healthy Smokers
if (length(subset(Sankey, Site == "Saliva" & Smoking == "Smokers" & Group == "Healthy")$Genus) <1 ) {
  Saliva = ""
}  else {
  Saliva = subset(Sankey, Site == "Saliva" & Smoking == "Smokers" & Group == "Healthy")$Genus
}

if (length(subset(Sankey, Site == "Plaque" & Smoking == "Smokers" & Group == "Healthy")$Genus) <1 ) {
  Plaque = ""
}  else {
  Plaque = subset(Sankey, Site == "Plaque" & Smoking == "Smokers" & Group == "Healthy")$Genus
}

if (length(subset(Sankey, Site == "Subgingival" & Smoking == "Smokers" & Group == "Healthy")$Genus) <1 ) {
  Subgingival = ""
}  else {
  Subgingival = subset(Sankey, Site == "Subgingival" & Smoking == "Smokers" & Group == "Healthy")$Genus
}
HS_Saliva <- Saliva
HS_Plaque <- Plaque
HS_Subgingival <- Subgingival
rm(Saliva, Plaque, Subgingival)


# Perio Non Smokers
if (length(subset(Sankey, Site == "Saliva" & Smoking == "Non Smokers" & Group == "Healthy")$Genus) <1 ) {
  Saliva = ""
}  else {
  Saliva = subset(Sankey, Site == "Saliva" & Smoking == "Non Smokers" & Group == "Healthy")$Genus
}

if (length(subset(Sankey, Site == "Plaque" & Smoking == "Non Smokers" & Group == "Healthy")$Genus) <1 ) {
  Plaque = ""
}  else {
  Plaque = subset(Sankey, Site == "Plaque" & Smoking == "Non Smokers" & Group == "Healthy")$Genus
}

if (length(subset(Sankey, Site == "Subgingival" & Smoking == "Non Smokers" & Group == "Healthy")$Genus) <1 ) {
  Subgingival = ""
}  else {
  Subgingival = subset(Sankey, Site == "Subgingival" & Smoking == "Non Smokers" & Group == "Healthy")$Genus
}

HNS_Saliva <- Saliva
HNS_Plaque <- Plaque
HNS_Subgingival <- Subgingival
rm(Saliva, Plaque, Subgingival)


# Chart
Plot_HNS <- venn.diagram(
                          x = list(HNS_Saliva, HNS_Plaque),
                          main = "Healthy - Non smokers",
                          category.names = c("Saliva", "Plaque"),
                          filename = NULL, 
                          output = TRUE ,
                          scaled = TRUE, 
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col=c("#440154ff", '#21908dff'),
                          #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                          cex = 0.5,
                          fontfamily = "sans",
                          cat.cex = 0.3,
                          cat.default.pos = "outer",
                          main.fontfamily = "Avenir",
                          cat.fontfamily = "Avenir")


# Chart
Plot_HS <- venn.diagram(
                          x = list(HS_Saliva, HS_Plaque),
                          main = "Healthy - Smokers",
                          category.names = c("Saliva", "Plaque"),
                          filename = NULL, 
                          output = TRUE ,
                          scaled = TRUE, 
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col=c("#440154ff", '#21908dff'),
                          #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                          cex = 0.5,
                          fontfamily = "sans",
                          cat.cex = 0.3,
                          cat.default.pos = "outer",
                          main.fontfamily = "Avenir",
                          cat.fontfamily = "Avenir")



Plot_PNS <- venn.diagram(
                          x = list(PNS_Saliva, PNS_Plaque, PNS_Subgingival),
                          main = "Periodontitis - Non smokers",
                          category.names = c("Saliva", "Plaque", "Subgingival"),
                          filename = NULL, 
                          output = TRUE ,
                          scaled = TRUE, 
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col=c("#440154ff", '#21908dff', "#9078ab"),
                          #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                          cex = 0.5,
                          fontfamily = "sans",
                          cat.cex = 0.3,
                          cat.default.pos = "outer",
                          main.fontfamily = "Avenir",
                          cat.fontfamily = "Avenir")





Plot_PS <- venn.diagram(
                          x = list( PS_Plaque, PS_Subgingival),
                          main = "Periodontitis - Smokers",
                          category.names = c( "Plaque", "Subgingival"),
                          filename = NULL, 
                          output = TRUE ,
                          scaled = TRUE, 
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col=c('#21908dff', "#9078ab"),
                          #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                          cex = 0.5,
                          fontfamily = "sans",
                          cat.cex = 0.3,
                          cat.default.pos = "outer",
                          main.fontfamily = "Avenir",
                          cat.fontfamily = "Avenir")






ggsave(Plot_HNS, 
       device = cairo_pdf,
       file="~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Healthy - Non smokers.pdf")


ggsave(Plot_HS, 
       device = cairo_pdf,
       file="~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Healthy - Smokers.pdf")

ggsave(Plot_PNS, 
       device = cairo_pdf,
       file="~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Perio - Non smokers.pdf")


ggsave(Plot_PS, 
       device = cairo_pdf,
       file="~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Perio - Smokers.pdf")





