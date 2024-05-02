library(phyloseq)
library(ggpubr)
library(ggplot2)
library(microbiome)
library(plyr)
library(reshape2)

# Load my 16S data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

Perio <- tax_glom(psm_perio, taxrank = "Phylum")
Healthy <- tax_glom(psm_healthy, taxrank = "Phylum")


MySite <- c("Saliva", "Plaque", "Subgingival")
res <- data.frame()

# Cycle for Periodontitis
for (i in MySite) {
  foo <- subset_samples(Perio, Site == i)
  foo <- transform(foo, transform = "compositional")
  foo <- psmelt(foo)
  foo$Sex <- factor(foo$Sex, levels = c("female", "male"), labels = c("F", "M"))
 foo <- subset(foo, Abundance >0.01)
  tmp <- compare_means(Abundance ~ Sex, group.by = "Phylum", method = "wilcox.test", p.adjust.method = "fdr", data = foo)
  x <- ddply(foo, c("Sex", "Phylum"), summarize,
             avg = mean(Abundance, na.rm = TRUE),
              N = sum(!is.na(Abundance)))
  x <- dcast(x, Phylum ~ Sex, value.var = "avg")
  tmp <- merge(tmp, x)
  tmp$Site <- i
  tmp$Periodontitis <- "Periodontitis"
  res <- rbind(res, tmp)
  rm(foo, tmp, x)
}

# Cycle for Healthy
for (i in MySite) {
  foo <- subset_samples(Healthy, Site == i)
  foo <- transform(foo, transform = "compositional")
  foo <- psmelt(foo)
  foo$Sex <- factor(foo$Sex, levels = c("female", "male"), labels = c("F", "M"))
  foo <- subset(foo, Abundance >0.01)
  tmp <- compare_means(Abundance ~ Sex, group.by = "Phylum", method = "wilcox.test", p.adjust.method = "fdr", data = foo)
  x <- ddply(foo, c("Sex", "Phylum"), summarize,
             avg = mean(Abundance, na.rm = TRUE),
             N = sum(!is.na(Abundance)))
  x <- dcast(x, Phylum ~ Sex, value.var = "avg")
  tmp <- merge(tmp, x)
  tmp$Site <- i
  tmp$Periodontitis <- "Healthy"
  res <- rbind(res, tmp)
  rm(foo, tmp, x)
}

df <- melt(res, id.vars = c("Phylum", "Site", "Periodontitis"), measure.vars = c("F", "M"))
names(df) <- c("Phylum", "Site", "Periodontitis", "Sex", "Abundance")

df$Site <- factor(df$Site, levels = c("Saliva", "Plaque", "Subgingival"), labels = c("S", "P", "sG"))


Healthy <- subset(df, Periodontitis == "Healthy")
Perio <- subset(df, Periodontitis == "Periodontitis")



MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 2/0.8,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(fill = NA, colour = NA))


p1 <- ggplot(Healthy, aes(x=Sex, y=Phylum, fill = Abundance))+
              geom_tile(color = "gray15")+
              facet_grid(~ Site)+
              scale_fill_distiller(palette = "BuGn", direction = 1)+
              MyTheme +
              labs(x="", y="", title = "Healthy")+
              theme(plot.margin = unit(c(0,0,-2,-2), "mm"),
                  panel.border = element_rect(colour = NA),
                  panel.spacing = unit(c(-2,-2), "mm"))


p2 <- ggplot(Perio, aes(x=Sex, y=Phylum, fill = Abundance))+
                geom_tile(color = "gray15")+
                facet_grid(~ Site)+
                scale_fill_distiller(palette = "BuGn", direction = 1)+
                MyTheme +
                labs(x="", y="", title = "Periodontitis")+
                theme(plot.margin = unit(c(0,0,-2,-2), "mm"),
                    panel.border = element_rect(colour = NA),
                    panel.spacing = unit(c(-2,-2), "mm"))

ggarrange(p1, p2)

ggsave(filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Heatmaps.eps",
        plot = last_plot(), width = 6, height = 5.5)
