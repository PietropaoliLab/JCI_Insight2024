library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiome)
# Heatmap.R
# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Phylum")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")


df <- psmelt(pseq_genus)
df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)

strata <- unique(df$Phylum)
res <- data.frame()
counter <- 1
for (i in strata) {
  
  Sex <- glm(AbundanceLog10 ~ Sex, data = subset(df, Phylum == i))
  Site <- glm(AbundanceLog10 ~ Site, data = subset(df, Phylum == i))
  Age <- glm(AbundanceLog10 ~ AgeYrs, data = subset(df, Phylum == i))
  Smoke <- glm(AbundanceLog10 ~ Smoking, data = subset(df, Phylum == i))
  Size <- glm(AbundanceLog10 ~ size_MB, data = subset(df, Phylum == i))
  Study <- glm(AbundanceLog10 ~ BioProject, data = subset(df, Phylum == i))
  #Race <- glm(AbundanceLog10 ~ Race, data = subset(df, Phylum == i))
  Country <- glm(AbundanceLog10 ~ Geo, data = subset(df, Phylum == i))
  
  
  
  # Extracting P values
  Psex <- coef(summary(Sex))[2,4]
  Psite <- coef(summary(Site))[2,4]
  Page <- coef(summary(Age))[2,4]
  Psmoke <- coef(summary(Smoke))[2,4]
  Psize <- coef(summary(Size))[2,4]
  Pstudy <- coef(summary(Study))[2,4]
  #Prace <- coef(summary(Race))[2,4]
  Pcountry <- coef(summary(Country))[2,4]
  
  foo <- data.frame("Phylum" = i, Psex, Psite, Page, Psmoke, Psize, Pstudy, "Prace", Pcountry, "Population" = "Healthy")
  res <- rbind(res, foo)
  rm(Psex, Psite, Page, Psmoke, Psize, Pstudy, Prace, Pcountry, Sex, Site, Smoke, Age, Size, Study, Race, Country)
  
  cat(paste0("Done: ", counter, " of ", length(strata), "         \r"))
  counter <- counter+1
}







# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Phylum")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")


df <- psmelt(pseq_genus)
df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)

strata <- unique(df$Phylum)


counter <- 1
for (i in strata) {
  
  Sex <- glm(AbundanceLog10 ~ Sex, data = subset(df, Phylum == i))
  Site <- glm(AbundanceLog10 ~ Site, data = subset(df, Phylum == i))
  Age <- glm(AbundanceLog10 ~ AgeYrs, data = subset(df, Phylum == i))
  Smoke <- glm(AbundanceLog10 ~ Smoking, data = subset(df, Phylum == i))
  Size <- glm(AbundanceLog10 ~ size_MB, data = subset(df, Phylum == i))
  Study <- glm(AbundanceLog10 ~ BioProject, data = subset(df, Phylum == i))
  #Race <- glm(AbundanceLog10 ~ Race, data = subset(df, Phylum == i))
  Country <- glm(AbundanceLog10 ~ Geo, data = subset(df, Phylum == i))
  
  
  
  # Extracting P values
  Psex <- coef(summary(Sex))[2,4]
  Psite <- coef(summary(Site))[2,4]
  Page <- coef(summary(Age))[2,4]
  Psmoke <- coef(summary(Smoke))[2,4]
  Psize <- coef(summary(Size))[2,4]
  Pstudy <- coef(summary(Study))[2,4]
  #Prace <- coef(summary(Race))[2,4]
  Pcountry <- coef(summary(Country))[2,4]
  
  foo <- data.frame("Phylum" = i, Psex, Psite, Page, Psmoke, Psize, Pstudy, "Prace", Pcountry, "Population" = "Periodontitis")
  res <- rbind(res, foo)
  rm(Psex, Psite, Page, Psmoke, Psize, Pstudy, Prace, Pcountry, Sex, Site, Smoke, Age, Size, Study, Race, Country)
  
  cat(paste0("Done: ", counter, " of ", length(strata), "         \r"))
  counter <- counter+1
}



healthy <- na.omit(subset(res, Population == "Healthy"))
rownames(healthy) <- healthy$Phylum
healthy <- healthy[,c(1:7, 9)]
colnames(healthy) <- c("Phylum","Sex", "Site", "Age", "Smoke", "Lib.Size", "Study", "Country")


perio <- na.omit(subset(res, Population == "Periodontitis"))
rownames(perio) <- perio$Phylum
perio <- perio[,c(1:7, 9)]
colnames(perio) <- c("Phylum","Sex", "Site", "Age", "Smoke", "Lib.Size", "Study", "Country")





library(reshape2)
helathy_melt <- reshape2::melt(healthy)

perio_melt <- reshape2::melt(perio)

library(ggplot2)
library(dplyr)

MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "left",
        legend.key = element_blank(),
        legend.background=element_blank(),
        plot.margin=unit(c(0.5,0,0.5,0), "cm"))


p <- ggplot(data = helathy_melt) + 
         geom_tile(data = helathy_melt %>% filter(value > 0.05), aes(x=variable, y=Phylum), fill="#f7f7f7", color = "Gray15")+
         geom_tile(data = helathy_melt %>% filter(value <= 0.05), aes(x=variable, y=Phylum, fill=value), color = "Gray15")+
  scale_fill_distiller(palette = "BuGn", direction = -1)+
      labs(title = "Healthy", x="", y="", fill = "P-value")+
          MyTheme+
        coord_fixed(ratio=2/1.5)+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

p

p1 <- ggplot(data = perio_melt) + 
  geom_tile(data = perio_melt %>% filter(value > 0.05), aes(x=variable, y=Phylum), fill="#f7f7f7", color = "Gray15")+
  geom_tile(data = perio_melt %>% filter(value <= 0.05), aes(x=variable, y=Phylum, fill=value), color = "Gray15")+
  scale_fill_distiller(palette = "BuGn", direction = -1)+
  labs(title = "Periodontitis", x="", y="", fill = "P-value")+
  MyTheme+
  coord_fixed(ratio=2/1.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5,-0.5), "cm"),
        legend.position = "none")

p1

library(gridExtra)
library(grid)
pdf(file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Heatmap Association.pdf")
grid.draw(cbind(ggplotGrob(p), ggplotGrob(p1), size="last"))
dev.off()

