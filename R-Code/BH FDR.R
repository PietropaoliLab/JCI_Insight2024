library(phyloseq) # also the basis of data object. Data analysis and visualisation

# Load my data
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# This method merges species that have the same taxonomy at a certain taxaonomic rank
pseq_genus <- phyloseq::tax_glom(psm_perio, taxrank = "Genus")

pseq_genus@sam_data$Site <- factor(pseq_genus@sam_data$Site)
pseq_genus@sam_data$Site <- relevel(pseq_genus@sam_data$Site, ref = "Saliva")

MyPhylum <- data.frame(tax_table(pseq_genus)[,"Genus"])
taxa_names(pseq_genus)
taxa_names(pseq_genus) <- MyPhylum$Genus
taxa_names(pseq_genus)

MySite <- unique(pseq_genus@sam_data$Site)

res <- data.frame()

for (x in MySite) {
  

tmp <- subset_samples(pseq_genus, Site == x)


df <- psmelt(tmp)
rm(tmp)

df <- abundances(df)
df$AbundanceLog10 <- log10(df$Abundance + 1)



strata <- unique(df$Genus)
for (i in strata) {
  p <- wilcox.test(AbundanceLog10 ~ Sex, data = subset(df, Genus == i))$p.value
  foo <- data.frame("Genus"=i, p, "Site"=x)
  res <- rbind(res, foo)
  cat(paste0(i, ": ", round(p,3), "\n"))
  rm(p, foo)
}
}

res$p_adjFDR <- p.adjust(res$p, method =  "BH")


library(ggplot2)
ggplot(res, aes(x=Genus, y=p_adjFDR))+
  geom_col()+
  facet_null(~Site)

attach(res)
newdata <- res[order(Site, p_adjFDR),]
detach(res)

newdata$ord <- as.factor(seq(1:nrow(newdata)))
newdata <- na.omit(newdata)

# plot everything
ggplot(newdata, aes(x = ord, y=p_adjFDR, fill=Site)) +  
  geom_bar(stat="identity")+
  geom_hline(yintercept = 0.05, size= .2, color = "black")+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.05, 0.8))


