library("microbial")
library(phyloseq)
library(microbiome)

library(gridExtra) 


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


res <- subset(res, p.value < 0.05 & LDAscore >2 & rank == "Genus")
res$y <- ifelse(res$direction == "male", -res$LDAscore, res$LDAscore) 
res$direction <- factor(res$direction, levels = c("female", "male"), labels = c("F", "M"))

res$tax <- gsub(".*_", "", res$tax)




MyTheme <- theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        #aspect.ratio = 1/2,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA))

MySite <- c("Saliva", "Plaque", "Subgingival")
MyPlots <- list() 
n <- 1

for (i in MySite) {
  
# Healthy
foo <- subset(res, Site ==i & Group == "Healthy")
p <- ggplot(foo, aes(x=reorder(tax, y), y=y, fill = direction, group = direction))+
          geom_hline(yintercept = 0, size = .3, linetype =1)+
          geom_bar(stat="identity", width = .6)+
          scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
          labs(title = "Healthy", subtitle = i)+
          coord_flip()+ 
          scale_y_continuous(limits=c(-6,6))+
          facet_grid(Group ~ Site, scales = "free", space = "free")+
          labs(y="LDA", x="")+
          MyTheme

  MyPlots[[n]] <- p # Plot Storing 
  n = n + 1 # Increment counter
  rm(foo, p) # Cleaning

  ##
# Periodontitis
foo <- subset(res, Site ==i & Group == "Periodontitis")
p <- ggplot(foo, aes(x=reorder(tax, y), y=y, fill = direction, group = direction))+
          geom_hline(yintercept = 0, size = .3, linetype =1)+
          geom_bar(stat="identity", width = .6)+
          scale_fill_manual(values = c("#F88B9D", "#8ECEFD"))+
          labs(title = "Periodontitis", subtitle = i)+
          coord_flip()+ 
          scale_y_continuous(limits=c(-6,6))+
          facet_grid(Group ~ Site, scales = "free", space = "free")+
          labs(y="LDA", x="")+
          MyTheme

  MyPlots[[n]] <- p # Plot Storing 
  n = n + 1 # Increment counter
  rm(foo, p) # Cleaning  
}

MyPlots <- MyPlots[1:5]

cowplot::plot_grid(plotlist=MyPlots,
                   align = "v",
                   nrow = 3,
                   rel_heights = c(.9, 1.5, .8))


