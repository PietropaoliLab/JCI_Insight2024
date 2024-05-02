library(phyloseq)
library(microbiome)
library(microbial)
library(extrafont)
library(showtext)
library(patchwork)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(ghibli)



loadfonts(device = "pdf")

# Load my 16S data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")


psm_perio <- tax_glom(psm_perio, taxrank = "Phylum")
psm_perio <- transform(psm_perio)
psm_perio@sam_data$Sex <- factor(psm_perio@sam_data$Sex, levels = c("male", "female"), labels = c("M", "F"))


psm_healthy <- tax_glom(psm_healthy, taxrank = "Phylum")
psm_healthy <- transform(psm_healthy)
psm_healthy@sam_data$Sex <- factor(psm_healthy@sam_data$Sex, levels = c("male", "female"), labels = c("M", "F"))




# My theme
MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Avenir"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(family = "Avenir", color = "black", size = 7),
        aspect.ratio = 1/15,
        legend.key = element_blank(),
        legend.background=element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Avenir", color = "black"))

# MyPalette
nb.cols <- 15
MyPalette <- colorRampPalette(c("#000000", "#2B4279", "#65428C", "#DFA6ED","#A1378B", "#D4237A",
                                "#F62B4C", "#FF5800", "#FFD286", "#E7E7E7"))(nb.cols)

# Setting my variables
Sites <- c("Saliva", "Plaque", "Subgingival")
MyPlots <- list()
res <- data.frame()

# # ================================
# #
# # Clycle for periodontitis
# 
# # ================================
# Within Sites
for (s in Sites) {
  # s = "Saliva"
  foo <- subset_samples(psm_perio, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  # Within dataset
  for (d in datasets){
    # d = "PRJNA321534"
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Print some info on screen
    cat(paste0("Analyzing: ", d, "-", s, "\n"))
    
    # Operate my analysis
    tmp <- plotbar(ps, level="Phylum", group="Sex", return = TRUE)
    
    tmp$BioProject <- d
    tmp$Site <- s
    tmp$Group <- "Periodontitis"
    
    # My index
    i <- paste0(d, " - ", "Periodontitis", " - ", s)
    
    # Compute average of Phyla abundance
    RelativeAbundancePlot <- ddply(tmp, c("Phylum", "Sex", "BioProject", "Site"), summarise,
                                   Abn = mean(Abundance))
    
    RelativeAbundancePlot$LongLabel <- i
    
    # Generating my ggplot
    p <- ggplot(RelativeAbundancePlot, aes(fill=Phylum, y=Abn, x=Sex)) + 
                    geom_bar(position="fill", stat="identity", width = 1)+
                       scale_fill_manual(values = MyPalette)+
                    #  scale_fill_brewer(palette = "Spectral", direction = 1)+
                              #facet_grid(Site~BioProject, drop = TRUE)+
                              labs(x="", y="Relative abundance (%)", fill ="Phylum", subtitle = i)+
                              coord_flip(expand = 0)+
                              MyTheme
    
    # Store my plot
    MyPlots[[i]] <- p
    
    res <- rbind(res, RelativeAbundancePlot)
    
    # Cleaning
    rm(tmp, p, RelativeAbundancePlot, i)
    
  }
}
rm(foo, ps, d, s)


#MyPlots
res$LongLabel <- factor(res$LongLabel, levels = c("PRJNA321534 - Periodontitis - Saliva", 
                                                  "PRJNA774981 - Periodontitis - Saliva",
                                                  "PRJNA774299 - Periodontitis - Saliva", 
                                                  "PRJNA321534 - Periodontitis - Plaque", 
                                                  "PRJNA324274 - Periodontitis - Plaque",
                                                  "PRJEB6047 - Periodontitis - Plaque", 
                                                  "PRJNA477241 - Periodontitis - Subgingival", 
                                                  "PRJNA773202 - Periodontitis - Subgingival",
                                                  "PRJEB6047 - Periodontitis - Subgingival"))



Perio <- ggplot(res, aes(x = Sex, y=Abn, fill = Phylum))+
              geom_bar(position="fill", stat="identity", width = 1)+
              facet_grid(LongLabel~.,)+
              coord_flip(expand = 0)+
              scale_fill_manual(values = MyPalette)+
              MyTheme+
              theme(panel.spacing.y = unit(0, "lines"),
                    strip.text.y = element_text(angle = 0),
                    strip.placement = "outside")


# ==============================================================================================
# 
#                     HEALTHY
# 
# ==============================================================================================

# MyPalette
# nb.cols <- 15
# MyPalette <- colorRampPalette(ghibli_palettes$MarnieMedium2)(nb.cols)
# MyPalette <- rev(MyPalette)
res <- data.frame()
# # ================================
# #
# # Clycle for healthy
# 
# # ================================
# Within Sites
for (s in Sites) {
  # s = "Saliva"
  foo <- subset_samples(psm_healthy, Site == s)
  
  # Get bioproject with specific site 
  datasets <- unique(foo@sam_data$BioProject)
  
  # Within dataset
  for (d in datasets){
    # d = "PRJNA321534"
    ps <- subset_samples(foo, BioProject == d) # Change here the dataframe
    
    # Print some info on screen
    cat(paste0("Analyzing: ", d, "-", s, "\n"))
    
    # Operate my analysis
    tmp <- plotbar(ps, level="Phylum", group="Sex", return = TRUE)
    
    tmp$BioProject <- d
    tmp$Site <- s
    tmp$Group <- "Healthy"
    
    # My index
    i <- paste0(d, " - ", "Healthy", " - ", s)
    
    # Compute average of Phyla abundance
    RelativeAbundancePlot <- ddply(tmp, c("Phylum", "Sex", "BioProject", "Site"), summarise,
                                   Abn = mean(Abundance))
    
    RelativeAbundancePlot$LongLabel <- i
    
    # Generating my ggplot
    p <- ggplot(RelativeAbundancePlot, aes(fill=Phylum, y=Abn, x=Sex)) + 
      geom_bar(position="fill", stat="identity", width = 1)+
      scale_fill_manual(values = MyPalette)+
      #  scale_fill_brewer(palette = "Spectral", direction = 1)+
      #facet_grid(Site~BioProject, drop = TRUE)+
      labs(x="", y="Relative abundance (%)", fill ="Phylum", subtitle = i)+
      coord_flip(expand = 0)+
      MyTheme
    
    # Store my plot
    MyPlots[[i]] <- p
    
    res <- rbind(res, RelativeAbundancePlot)
    
    # Cleaning
    rm(tmp, p, RelativeAbundancePlot, i)
    
  }
}
rm(foo, ps, d, s)


#MyPlots

dput(unique(res$LongLabel))
res$LongLabel <- factor(res$LongLabel, levels = c( "PRJNA774981 - Healthy - Saliva", 
                                                   "PRJNA774299 - Healthy - Saliva", 
                                                   "PRJNA321534 - Healthy - Saliva",
                                                   "PRJNA321534 - Healthy - Plaque", 
                                                   "PRJEB6047 - Healthy - Plaque", 
                                                   "PRJNA773202 - Healthy - Subgingival",
                                                   "PRJEB6047 - Healthy - Subgingival" ))



Healthy <- ggplot(res, aes(x = Sex, y=Abn, fill = Phylum))+
                      geom_bar(position="fill", stat="identity", width = 1)+
                      facet_grid(LongLabel~.,)+
                      coord_flip(expand = 0)+
                      scale_fill_manual(values = MyPalette)+
                      MyTheme+
                      theme(panel.spacing.y = unit(0, "lines"),
                            strip.text.y = element_text(angle = 0),
                            strip.placement = "outside")

Healthy

library(patchwork)
RelativeAbundance <- Healthy/Perio + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
RelativeAbundance


ggsave(plot = RelativeAbundance, 
       device = cairo_pdf,
        width = 9, 
       # height = 11, 
       units = "in",
       filename = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Images/Relative Abundance Figure 1.pdf")





