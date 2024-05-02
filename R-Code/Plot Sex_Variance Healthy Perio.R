library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiome)

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")


df <- psmelt(pseq_genus)
df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)

strata <- unique(df$Genus)
res <- data.frame()
counter <- 1
for (i in strata) {
  
  fit <- glm(AbundanceLog10 ~ Sex + Site, data = subset(df, Genus == i))

  # Abundance mean
  Abn <- fit[["coefficients"]][["(Intercept)"]]
  
  
  # Compute McFadden's R-squared for model which ranges from 0 to just under 1, 
  # with higher values indicating a better model fit.
  # R-squared represents the proportion of the variance in the response variable
  # that can be explained by the predictor variables in a regression model.
  R <- with(summary(fit), 1 - deviance/null.deviance)
 
  
  
  # Extracting P values
  Psex <- coef(summary(fit))[2,4]
  
  foo <- data.frame("Species" = i, "Variance"= R, "P_sex" = Psex, "Abn" = Abn, "Population" = "Healthy")
  res <- rbind(res, foo)
  rm(foo, fit, R, Psex)
  
  cat(paste0("Done: ", counter, " of ", length(strata), "         \r"))
  counter <- counter+1
}



# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")


df <- psmelt(pseq_genus)
df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)

strata <- unique(df$Genus)
counter <- 1

for (i in strata) {
  
  fit <- glm(AbundanceLog10 ~ Sex + Site, data = subset(df, Genus == i))
  
  # Abundance mean
  Abn <- fit[["coefficients"]][["(Intercept)"]]
  
  
  # Compute McFadden's R-squared for model which ranges from 0 to just under 1, 
  # with higher values indicating a better model fit.
  # R-squared represents the proportion of the variance in the response variable
  # that can be explained by the predictor variables in a regression model.
  R <- with(summary(fit), 1 - deviance/null.deviance)
  
  
  
  # Extracting P values
  Psex <- coef(summary(fit))[2,4]
  
  foo <- data.frame("Species" = i, "Variance"= R, "P_sex" = Psex, "Abn" = Abn, "Population" = "Periodontitis")
  res <- rbind(res, foo)
  rm(foo, fit, R, Psex)
  
  cat(paste0("Done: ", counter, " of ", length(strata), "         \r"))
  counter <- counter+1
}





# Compute FDR
res$FDR_sex <- p.adjust(res$P_sex, method = "BH")





res$FDR <- ifelse(res$FDR_sex < .05, "FDR<.05", "NS")
res$MyLabs <- ifelse(res$FDR == "FDR<.05", res$Species, NA)

library(ggplot2)
library("ggrepel")
library(plyr)
res <- res[order(res$FDR, decreasing=TRUE), ]

MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 2/1,
        legend.key = element_blank(),
        legend.background=element_blank())

ggplot(subset(res, !is.na(FDR)), aes(x=Variance, y=Abn, group = FDR))+
  geom_point(aes(color = FDR, shape = Population), size = 3.5)+
  scale_color_manual(values = c("#ca0020", "#d9d9d9"))+
  scale_shape_manual(values=c(16, 17))+
  #coord_fixed(xlim = c(0, .45), ylim = c(0, .45))+
  geom_label_repel(aes(label = MyLabs), size = 3, fill="white", seed = 101,
                   box.padding = unit(1.5, "lines"),
                   point.padding = unit(0.3, "lines"))+
  labs(x="Variance explained by sex", y=expression(paste(Log[10]," Abundance", sep = "")))+
  MyTheme+
  theme(legend.position = "right",
        legend.direction = "vertical")









