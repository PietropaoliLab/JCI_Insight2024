library(phyloseq)
library(microbiome) # also the basis of data object. Data analysis and visualisation
library(svMisc)

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM PlyloSeq Object for analysys.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_ps, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")

MyGenus <- data.frame(tax_table(pseq_genus)[,"Genus"])
taxa_names(pseq_genus)


df <- psmelt(pseq_genus)
df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)

strata <- unique(df$Genus)
MySite <- unique(as.character(df$Site))
GenusCounter <-  1 

res <- data.frame()

  for (i in strata) {
    for (site in MySite) {
  fitSex <- glm(AbundanceLog10 ~ Sex, data = subset(df, Genus ==i & Site == site))
  #fitSite <- glm(AbundanceLog10 ~ Site, data = subset(df, Genus ==i & Site == site))
  fitStudy <- glm(AbundanceLog10 ~ BioProject, data = subset(df, Genus ==i & Site == site))
  fitPerio <- glm(AbundanceLog10 ~ Periodontitis, data = subset(df, Genus ==i & Site == site))
  
  
  # Abundance mean
  AbnSex <- fitSex[["coefficients"]][["(Intercept)"]]
  #AbnSite <- fitSite[["coefficients"]][["(Intercept)"]]
  AbnStudy <- fitStudy[["coefficients"]][["(Intercept)"]]
  AbnPerio <- fitPerio[["coefficients"]][["(Intercept)"]]
  
  
  
  # Compute McFadden's R-squared for model which ranges from 0 to just under 1, 
  # with higher values indicating a better model fit.
  # R-squared represents the proportion of the variance in the response variable
  # that can be explained by the predictor variables in a regression model.
  Rsex <- with(summary(fitSex), 1 - deviance/null.deviance)
  #Rsite <- with(summary(fitSite), 1 - deviance/null.deviance)
  Rstudy <- with(summary(fitStudy), 1 - deviance/null.deviance)
  Rperio <- with(summary(fitPerio), 1 - deviance/null.deviance)
  
  
  foo <- data.frame("VarianceSex"= Rsex, 
                    #"VarianceSite"= Rsite, 
                    "VarianceStudy"= Rstudy, "VariancePerio"= Rperio,
                    "Species" = i,
                    "Site" = site, 
                    AbnSex, 
                    #AbnSite, 
                    AbnStudy, AbnPerio)
  res <- rbind(res, foo)
  rm(foo, fitSex, fitSite, fitStudy, fitPerio, Rsex, Rsite, Rstudy,Rperio, AbnSex, AbnSite, AbnStudy, AbnPerio)
  
  cat("\r")
  cat(paste0("Genera:", GenusCounter, " of ",length(strata) * length(MySite), " - Site: ", site, "                     \r") )
  GenusCounter <- GenusCounter + 1 
  
} 
  } 



# Generating trimmed mean 10%
f <- function(x) mean(x, trim=0.1)
res$Abundance <- apply(res[6:8], 1, f)

library(ggplot2)
MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "top")

ggplot(res, aes(x=VarianceStudy, y=VariancePerio))+
  geom_point(aes(size = Abundance, fill = Site), shape = 21, color = "black", alpha=.5)+
  #coord_fixed(xlim = c(0, .45), ylim = c(0, .45))+
  labs(x="Variance explained by study", y="Variance explained by periodontitis")+
  MyTheme


ggplot(res, aes(x=VarianceStudy, y=VarianceSite))+
  geom_point(aes(size = Abundance), shape = 21, fill = "#377eb8", color = "black", alpha=.3)+
  coord_fixed(xlim = c(0, .7), ylim = c(0, .7))+
  labs(x="Variance explained by study", y="Variance explained by sampling site")+
  MyTheme



GLMfit <- glm(AbundanceLog10 ~ Sex, data = subset(df, Genus =="Porphyromonas" & Site == "Saliva"))
LMfit <-   lm(AbundanceLog10 ~ Sex, data = subset(df, Genus =="Porphyromonas" & Site == "Saliva"))

Rsex <- with(summary(GLMfit), 1 - deviance/null.deviance)
Rsex
summary(LMfit)

