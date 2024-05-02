library(UpSetR)
library(xlsx)
library(extrafont)
library(showtext)
loadfonts(device = "pdf")


Combi <- read.xlsx(file = "AQ Projects/Microbiome/Gender features of periodontal microbiome/Sankey data.xlsx", sheetIndex = 1)

Combi$Agreement <- TRUE

Combi$Periodontitis <- ifelse(Combi$Group == "Periodontitis", 1, 0)
Combi$Healthy <- ifelse(Combi$Group == "Healthy", 1, 0)

Combi$Smokers <- ifelse(Combi$Smoking == "Smokers", 1, 0)
Combi$`Non Smokers` <- ifelse(Combi$Smoking == "Non Smokers", 1, 0)

Combi$Saliva <- ifelse(Combi$Site == "Saliva", 1, 0)
Combi$Plaque <- ifelse(Combi$Site == "Plaque", 1, 0)
Combi$Subgingival <- ifelse(Combi$Site == "Subgingival", 1, 0)

Combi <- Combi[, c(2,6:13)]



upsetPlot <- upset(Combi,
              nsets = 7,
              point.size = 3.5, 
              line.size = 1, 
              mainbar.y.label = "N. Genera identified as\ndifferent between M vs F", sets.x.label = "Genera per stratum")

upsetPlot

pdf(file="~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Upset.eps") # or other device
upsetPlot
dev.off()

p


ggsave('~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Upset.eps', 
       ggplotify::as.ggplot(upsetPlot),
       device = cairo_pdf,
       width = 8, 
       units = "in")

