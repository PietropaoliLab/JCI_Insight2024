load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM PlyloSeq Object for analysys.RData")

library(MicrobiotaProcess)
library(patchwork)
library(ggplot2)
library(plyr)

#set.seed(1024)
#ps <- rarefy_even_depth(psm_ps)
#ps
# rareres <- get_rarecurve(obj = ps, chunks=400)
load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/Rarecurve OBJ for rarefaction Plot.RData")

df <- as.data.frame(rareres$data)

cdata <- ddply(subset(df, Alpha == "Shannon"), c("readsNums","Sex", "Site", "Smoking", "Periodontitis"), summarise,
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


NSmoking <- ggplot(subset(cdata, Smoking == "No"), aes(readsNums, mean, group = Sex, color = Sex, fill = Sex))+
 # geom_vline(xintercept = 500, linetype = 2, size =.2)+
                  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = .2, colour = NA)+
                  geom_line(size=.75)+
                    facet_grid(Periodontitis~Site)+
                    labs(title = "Non smokers", x="Numbers of reads", y="Shannon index")+
                    scale_fill_manual(values=c("#8ECEFD", "#F88B9D"))+
                    scale_color_manual(values=c("#8ECEFD", "#F88B9D"))+
                  MyTheme +
                  coord_fixed(ylim = c(2.2, 5))+
                  theme(legend.position = c(.9,.15))


Smoking <- ggplot(subset(cdata, Smoking == "Yes"), aes(readsNums, mean, group = Sex, color = Sex, fill = Sex))+
  #geom_vline(xintercept = 500, linetype = 2, size =.2)+
                  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = .2, colour = NA)+
                  geom_line(size=.75)+
                  facet_grid(Periodontitis~Site)+
                    labs(title = "Smokers", x="Numbers of reads", y="Shannon index")+
                    scale_fill_manual(values=c("#8ECEFD", "#F88B9D"))+
                    scale_color_manual(values=c("#8ECEFD", "#F88B9D"))+
                    MyTheme +
                    coord_fixed(ylim = c(2.2, 5))+
                    theme(legend.position = c(.9,.15))


library(gridExtra)
grid.arrange(NSmoking, Smoking, nrow = 1) # Save as pdf 6.56*11.56

