library(phyloseq)
library(microbiome)
library(microbial)
library(extrafont)
library(showtext)
library(patchwork)
loadfonts(device = "pdf")

# preRef(ref_db = "silva", path=".")


load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

ps <- psm_perio
ps <- tax_glom(ps, taxrank = "Phylum")
#default normalize method is relative
phy <- transform(ps, transform = "compositional")

MySite <- unique(phy@sam_data$Site)

PerioRelativeAbn <- data.frame()

for (i in MySite) {
  
  foo <- subset_samples(phy, Site == i)
  tmp <- plotbar(foo, level="Phylum", group="Sex", return = TRUE)
  tmp$Site <- i
  tmp$Group <- "Periodontitis"
  PerioRelativeAbn <- rbind(PerioRelativeAbn, tmp)
  rm(foo, tmp)
}


# Healthy
ps <- psm_healthy
ps <- tax_glom(ps, taxrank = "Phylum")
#default normalize method is relative
phy <- transform(ps, transform = "compositional")

MySite <- unique(phy@sam_data$Site)

HealthyRelativeAbn <- data.frame()

for (i in MySite) {
  
  foo <- subset_samples(phy, Site == i)
  tmp <- plotbar(foo, level="Phylum", group="Sex", return = TRUE)
  tmp$Site <- i
  tmp$Group <- "Healthy"
  HealthyRelativeAbn <- rbind(HealthyRelativeAbn, tmp)
  rm(foo, tmp)
}


foo <- rbind(HealthyRelativeAbn, PerioRelativeAbn)
foo$Sex <- factor(foo$Sex, levels = c("female", "male"), labels = c("F", "M"))
foo$Site <- factor(foo$Site, levels = c("Saliva", "Plaque", "Subgingival"))

library(plyr)
RelativeAbundancePlot <- ddply(foo, c("Phylum", "Sex", "Site", "Group"), summarise,
                               Abn = mean(Abundance))

# x <- arrange(RelativeAbundancePlot,desc(Abn))
# MyTop10Phylum <- unique(x$Phylum)[1:10]
# 
# RelativeAbundancePlot <- subset(RelativeAbundancePlot, Phylum %in% MyTop10Phylum)


library(ggplot2)
MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Avenir"),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(family = "Avenir", color = "black", size = 7),
        aspect.ratio = 2/.6,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir", color = "black"))

# MyPalette <- c("#000000", "#2B4279", "#65428C", "#DFA6ED","#A1378B", "#D4237A",
#   "#F62B4C", "#FF5800", "#FFD286", "#E7E7E7" )

# MyPalette
nb.cols <- 15
MyPalette <- colorRampPalette(c("#000000", "#2B4279", "#65428C", "#DFA6ED","#A1378B", "#D4237A",
                                "#F62B4C", "#FF5800", "#FFD286", "#E7E7E7"))(nb.cols)


P_Phylum <- ggplot(RelativeAbundancePlot, aes(fill=Phylum, y=Abn, x=Sex)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = MyPalette)+
  #scale_fill_brewer(palette = "Spectral", direction = 1)+
  facet_grid(Group~Site)+
  labs(x="", y="Relative abundance (%)", fill ="Phylum")+
  coord_fixed(expand = 0)+
  MyTheme

P_Phylum

ggsave(plot = P_Phylum, 
       device = cairo_pdf,
       # width = 5.5, 
       # height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Relative_abundances_Phylum.pdf")





ps <- psm_perio
ps <- tax_glom(ps, taxrank = "Genus")
#default normalize method is relative
phy <- transform(ps, transform = "compositional")

MySite <- unique(phy@sam_data$Site)

PerioRelativeAbn <- data.frame()

for (i in MySite) {
  
  foo <- subset_samples(phy, Site == i)
  tmp <- plotbar(foo, level="Genus", group="Sex", return = TRUE)
  tmp$Site <- i
  tmp$Group <- "Periodontitis"
  PerioRelativeAbn <- rbind(PerioRelativeAbn, tmp)
  rm(foo, tmp)
}


# Healthy
ps <- psm_healthy
ps <- tax_glom(ps, taxrank = "Genus")
#default normalize method is relative
phy <- transform(ps, transform = "compositional")

MySite <- unique(phy@sam_data$Site)

HealthyRelativeAbn <- data.frame()

for (i in MySite) {
  
  foo <- subset_samples(phy, Site == i)
  tmp <- plotbar(foo, level="Genus", group="Sex", return = TRUE)
  tmp$Site <- i
  tmp$Group <- "Healthy"
  HealthyRelativeAbn <- rbind(HealthyRelativeAbn, tmp)
  rm(foo, tmp)
}


foo <- rbind(HealthyRelativeAbn, PerioRelativeAbn)
foo$Sex <- factor(foo$Sex, levels = c("female", "male"), labels = c("F", "M"))
foo$Site <- factor(foo$Site, levels = c("Saliva", "Plaque", "Subgingival"))

library(plyr)
RelativeAbundancePlotGenus <- ddply(foo, c("Genus", "Sex", "Site", "Group"), summarise,
                               Abn = mean(Abundance))

x <- arrange(RelativeAbundancePlotGenus,desc(Abn))
MyTop10Genus <- unique(x$Genus)[1:10]

RelativeAbundancePlotGenus <- subset(RelativeAbundancePlotGenus, Genus %in% MyTop10Genus)



P_Genus <- ggplot(RelativeAbundancePlotGenus, aes(fill=Genus, y=Abn, x=Sex)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "Set3")+
  facet_grid(Group~Site)+
  labs(x="", y="Relative abundance (%)", fill ="Genus")+
  MyTheme


ggsave(plot = P_Genus, 
       device = cairo_pdf,
       # width = 5.5, 
       # height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Relative_abundances_Genus.pdf")



P_Phylum+P_Genus

ggsave(plot = last_plot(), 
       device = cairo_pdf,
       # width = 7, 
       # height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Relative_abundances_Phylum_Genus.pdf")











# Beta Diversity
healthy<-normalize(psm_healthy)
p.h <- plotbeta(healthy, group="Site", distance = "bray", ellipse = FALSE, color = c("#7fc97f", "#1f78b4", "#d01c8b"))+
  MyTheme + labs(title = "Healthy", color = "Site") + theme(aspect.ratio = 1)
p.h$layers[[1]]$aes_params$alpha <- .3

perio<-normalize(psm_perio)
p.p <- plotbeta(perio, group="Site", distance = "bray", ellipse = FALSE, color = c("#7fc97f", "#1f78b4", "#d01c8b"))+
  MyTheme + labs(title = "Periodontitis", color = "Site") + theme(aspect.ratio = 1)
p.p$layers[[1]]$aes_params$alpha <- .3


pdf(file = "Microbiome/Gender features of periodontal microbiome/Results/Figures/Beta Diversity BRAY.pdf",
    width = 6, height = 3)
ggpubr::ggarrange(p.h, p.p, common.legend = TRUE)
dev.off()




load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

res <- data.frame()
MySite <- c("Saliva", "Plaque", "Subgingival")
for (i in MySite) {
  
  foo <- subset_samples(psm_healthy, Site ==i)
  tmp <- ldamarker(foo, group="Sex")
  tmp$Site = i
  tmp$Group = "Healthy"
  res <- rbind(res, tmp)
  rm(foo, tmp)
  
  foo <- subset_samples(psm_perio, Site ==i)
  tmp <- ldamarker(foo, group="Sex")
  tmp$Site = i
  tmp$Group = "Periodontitis"
  res <- rbind(res, tmp)
  rm(foo, tmp)
  
}


res <- subset(res, p.value < 0.05 & LDAscore >3)
res$y <- ifelse(res$direction == "male", -res$LDAscore, res$LDAscore) 
res$direction <- factor(res$direction, levels = c("male", "female"), labels = c("M", "F"))


ggplot(res, aes(x=reorder(tax, y), y=y, fill = direction, group = direction))+
  geom_hline(yintercept = 0, size = .3, linetype =1)+
  geom_bar(stat="identity")+
  coord_flip()+ 
  scale_y_continuous(limits=c(-6,6))+
  facet_grid(Group ~ Site, scales = "free", space = "free")+
  labs(y="LDA", x="")+
theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        
        legend.position = "right",
        #aspect.ratio = 1/2,
        legend.key = element_blank(),
        legend.background=element_blank())









