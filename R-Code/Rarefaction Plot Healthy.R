library(MicrobiotaProcess)
library(patchwork)
library(ggplot2)
library(plyr)
library(phyloseq)

# load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")

set.seed(1024)
# ps <- rarefy_even_depth(psm_healthy)
# ps
# HealthyRarefactionPlot <- get_rarecurve(obj = ps, chunks=400)

# Saving Rarefaction Plot Object for further use
# save(HealthyRarefactionPlot , file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Healthy Rarefaction Plot OBJ.Rdata")

# Loading Rarefaction Plot obj
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Healthy Rarefaction Plot OBJ.Rdata")

df <- as.data.frame(HealthyRarefactionPlot$data)

cdata <- ddply(subset(df, Alpha == "Shannon"), c("readsNums","Sex", "Site"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N))
cdata$Sex <- factor(cdata$Sex, levels = c("male", "female"), labels = c("Male", "Female"))

MyTheme <- theme_bw(12)+
  theme(axis.text=element_text(size=10, color = "black"), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey95"),
        aspect.ratio = 1,
        legend.title = element_blank())


p <- ggplot(cdata, aes(readsNums, mean, group = Sex, color = Sex, fill = Sex))+
  # geom_vline(xintercept = 500, linetype = 2, size =.2)+
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = .2, colour = NA)+
  geom_line(size=.75)+
  facet_grid(~Site)+
  labs(title = "Healthy individuals", x="Numbers of reads", y="Shannon index")+
  scale_fill_manual(values=c("#8ECEFD", "#F88B9D"))+
  scale_color_manual(values=c("#8ECEFD", "#F88B9D"))+
  MyTheme +
  coord_fixed(ylim = c(2.2, 5))+
  theme(legend.position = c(.9,.15))


pdf(file = "~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Results/Figures/Rarefaction Plot Healthy.pdf",
    width = 8.3, height = 6.25)
p

# Closing the graphical device
dev.off() 
p
