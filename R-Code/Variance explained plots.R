library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiome)

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_healthy, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")

MyGenus <- data.frame(tax_table(pseq_genus)[,"Genus"])
taxa_names(pseq_genus)
taxa_names(pseq_genus) <- MyGenus$Genus
taxa_names(pseq_genus)

df <- psmelt(pseq_genus)
df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)

strata <- unique(df$Genus)
res <- data.frame()

for (i in strata) {
  
  fitSex <- glm(AbundanceLog10 ~ Sex, data = subset(df, Genus ==i))
  fitSite <- glm(AbundanceLog10 ~ Site, data = subset(df, Genus ==i))
  fitStudy <- glm(AbundanceLog10 ~ BioProject, data = subset(df, Genus ==i))
  
  # Abundance mean
  AbnSex <- fitSex[["coefficients"]][["(Intercept)"]]
  AbnSite <- fitSite[["coefficients"]][["(Intercept)"]]
  AbnStudy <- fitStudy[["coefficients"]][["(Intercept)"]]
  
  
  # Compute McFadden's R-squared for model which ranges from 0 to just under 1, 
  # with higher values indicating a better model fit.
  # R-squared represents the proportion of the variance in the response variable
  # that can be explained by the predictor variables in a regression model.
  Rsex <- with(summary(fitSex), 1 - deviance/null.deviance)
  Rsite <- with(summary(fitSite), 1 - deviance/null.deviance)
  Rstudy <- with(summary(fitStudy), 1 - deviance/null.deviance)
  
  
  # Extracting P values
  Psex <- coef(summary(fitSex))[2,4]
  Psite <- coef(summary(fitSite))[2,4]
  Pstudy <- coef(summary(fitStudy))[2,4]
  
  
  
  foo <- data.frame("VarianceSex"= Rsex, "VarianceSite"= Rsite, "VarianceStudy"= Rstudy, "Species" = i, 
                    AbnSex, AbnSite, AbnStudy, "P_sex" = Psex, "P_site" = Psite, "P_study" = Pstudy)
  res <- rbind(res, foo)
  rm(foo, fitSex, fitSite, fitStudy, Rsex, Rsite, Rstudy, AbnSex, AbnSite, AbnStudy)
}

# Compute FDR
res$FDR_sex <- p.adjust(res$P_sex, method = "BH")
res$FDR_site <- p.adjust(res$P_site, method = "BH")
res$FDR_study <- p.adjust(res$P_study, method = "BH")

# Generating trimmed mean 10%
f <- function(x) mean(x, trim=0.1)
res$Abundance <- apply(res[5:7], 1, f)


res$FDR_site_study <- ifelse(res$P_site < 5e-5 & res$P_study < 5e-5, "FDR<5e-5", "NS")
res$FDR_site_sex <- ifelse(res$P_site < 5e-5 & res$P_sex < 5e-5, "FDR<5e-5", "NS")
res$FDR_study_sex <- ifelse(res$P_study < 5e-5 & res$P_sex < 5e-5, "FDR<5e-5", "NS")


library(ggplot2)

MyTheme <-   theme_bw(12)+
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        aspect.ratio = 1)


ggplot(subset(res, !is.na(FDR_site_sex)), aes(x=VarianceSex, y=VarianceSite))+
  geom_point(aes(color = FDR_site_sex), size = 3)+
  scale_color_manual(values = c("#d9d9d9"))+
  coord_fixed(xlim = c(0, 0.2), ylim = c(0, 0.6))+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(x="Variance explained by sex", y="Variance explained by sampling site")+
  MyTheme


ggplot(subset(res, !is.na(FDR_study_sex)), aes(x=VarianceSex, y=VarianceStudy))+
  geom_point(aes(color = FDR_study_sex), size = 3)+
  scale_color_manual(values = c("#d9d9d9"))+
  coord_fixed(xlim = c(0, 0.2), ylim = c(0, 0.8))+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  labs(x="Variance explained by sex", y="Variance explained by study")+
  MyTheme

  
  
  
  
  
