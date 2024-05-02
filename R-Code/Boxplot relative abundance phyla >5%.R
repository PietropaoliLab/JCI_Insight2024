library(phyloseq)
library(ggpubr)
library(ggplot2)
library(microbiome)
# Load my 16S data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

Perio <- tax_glom(psm_perio, taxrank = "Phylum")
Healthy <- tax_glom(psm_healthy, taxrank = "Phylum")

Perio <- transform(Perio, transform = "compositional")
Healthy <- transform(Healthy, transform = "compositional")

Perio <- psmelt(Perio)
Healthy <- psmelt(Healthy)

Perio <- subset(Perio, Abundance >0.05)
Healthy <- subset(Healthy, Abundance >0.05)

Perio$Sex <- factor(Perio$Sex, levels = c("female", "male"), labels = c("F", "M"))
Healthy$Sex <- factor(Healthy$Sex, levels = c("female", "male"), labels = c("F", "M"))


MyTheme <-   theme_bw(10)+
  theme(axis.text = element_text(colour = "black", size = 8),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 2/1.2,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.text.x = element_text(size = 7),
        strip.background = element_rect(fill = NA, colour = NA))

p1 <- ggboxplot(Perio, x="Sex", y="Abundance",
          color = "Sex", palette = c("#F88B9D", "#8ECEFD"),
          facet.by = "Phylum",
          add = "jitter", add.params = list(size = 1.5, alpha = .3, jitter =.1))+
          ylim(0,1)+
          labs(x="", y="Phylum relative abundance (%)", title = "Periodontitis")+
          stat_compare_means(aes(color = "Sex"), method = "wilcox.test", size = 2.5,
                             label = "p.format", label.y = .9, label.x = 1.25, show.legend = FALSE)+
          MyTheme
p1 <- facet(p1, facet.by = "Phylum", nrow = 2)



p2 <- ggboxplot(Healthy, x="Sex", y="Abundance",
                color = "Sex", palette = c("#F88B9D", "#8ECEFD"),
                facet.by = "Phylum",
                add = "jitter", add.params = list(size = 1.5, alpha = .3, jitter =.1))+
                ylim(0,1)+
                labs(x="", y="Phylum relative abundance (%)", title = "Healthy")+
                stat_compare_means(aes(color = "Sex"), method = "wilcox.test", size = 2.5,
                     label = "p.format", label.y = .9, label.x = 1.25, show.legend = FALSE)+
                 MyTheme
p2 <- facet(p2, facet.by = "Phylum", nrow = 2)

library(grid)
pdf(file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/BoxPlot Relative Abundance.pdf",
    width = 6, height = 8, onefile = TRUE)

grid.draw(rbind(ggplotGrob(p2), ggplotGrob(p1)))
dev.off()





