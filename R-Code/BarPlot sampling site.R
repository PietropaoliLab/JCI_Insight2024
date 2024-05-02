library(plyr)
library(ggplot2)
library(phyloseq)
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# Melting Phyloseq Objects
healthy <- psmelt(psm_healthy)
perio <- psmelt(psm_perio)

CombHealthy <- ddply(healthy, c("Site"), summarise,
               N = sum(!is.na(Site)))
CombHealthy$Perc <-  CombHealthy$N/sum(CombHealthy$N)
CombHealthy$Group <- "H"

CombPerio <- ddply(perio, c("Site"), summarise,
                     N = sum(!is.na(Site)))
CombPerio$Perc <-  CombPerio$N/sum(CombPerio$N)
CombPerio$Group <- "P"

DataPlot <- rbind(CombHealthy, CombPerio)

MyTheme <-   theme_minimal(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.key = element_blank(),
        aspect.ratio = 3/1,
        legend.background=element_blank() )

ggplot(DataPlot, aes(x= Group, y = Perc, fill = Site))+
         geom_col(position = "stack", alpha = .8, width = .8)+
  labs(x="", y="Sampling site (%)", fill ="")+
  scale_fill_manual(values = c("#7fc97f", "#1f78b4", "#d01c8b"))+
  MyTheme
       
       
       