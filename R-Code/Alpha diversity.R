# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling  

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM PlyloSeq Object for analysys.RData")


print(psm_ps)

# As is evident there is a large difference in the number of reads. Minimum is 7477 and maximum is 346214!! 
summary(sample_sums(psm_ps))


set.seed(101)  # This will help in reproducing the filtering and nomalisation. 
ps0.rar <- rarefy_even_depth(psm_ps, sample.size = 7477)

print(ps0.rar)
p.rar <- plot_taxa_prevalence(ps0.rar, "Family")
p.rar


# get the metadata out as seprate object
hmp.meta <- meta(ps0.rar)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)

# Add the rownames to diversity table
# Let us calculate diversity.
hmp.div <- alpha(psm_ps, index = "all")
hmp.div$Run <- rownames(hmp.div)

# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "Run")

# check the tables
colnames(div.df)

# Now use this data frame to plot 
p <- ggboxplot(subset(div.df, Smoking == "No"), 
               x = "Sex", 
               y = "diversity_shannon",
               fill = "Sex", 
               palette = c("#8ECEFD", "#F88B9D"),
               facet.by = c("Periodontitis","Site"))+
  stat_compare_means(label.y = 6, label = "p.format")+
  ggtitle("Non smoking")+
  labs(y="Shannon diversity")+
  scale_x_discrete(labels=c("male" = "M", "female" = "F"))+
  theme(legend.position = "none", aspect.ratio = 3.5)
print(p)


# Now use this data frame to plot 
p1 <- ggboxplot(subset(div.df, Smoking == "Yes"), 
               x = "Sex", 
               y = "diversity_shannon",
               fill = "Sex", 
               palette = c("#8ECEFD", "#F88B9D"),
               facet.by = c("Periodontitis","Site"))+
  stat_compare_means(label.y = 6, label = "p.format")+
  ggtitle("Smoking")+ 
  labs(y="Shannon diversity")+
  scale_x_discrete(labels=c("male" = "M", "female" = "F"))+
  theme(legend.position = "none", aspect.ratio = 3.5)
p1

library(gridExtra)
grid.arrange(p, p1, nrow=1)



# Now use this data frame to plot 
p2 <- ggboxplot(subset(div.df, Smoking == "No"), 
               x = "Sex", 
               y = "chao1",
               fill = "Sex", 
               palette = c("#8ECEFD", "#F88B9D"),
               facet.by = c("Periodontitis","Site"))+
  stat_compare_means(label.y = 650, label = "p.format")+
  ggtitle("Non smoking")+
  labs(y="Chao index")+
  scale_x_discrete(labels=c("male" = "M", "female" = "F"))+
  theme(legend.position = "none", aspect.ratio = 3.5)
print(p2)


# Now use this data frame to plot 
p3 <- ggboxplot(subset(div.df, Smoking == "Yes"), 
                x = "Sex", 
                y = "chao1",
                fill = "Sex", 
                palette = c("#8ECEFD", "#F88B9D"),
                facet.by = c("Periodontitis","Site"))+
  stat_compare_means(label.y = 650, label = "p.format")+
  ggtitle("Smoking")+ 
  labs(y="Chao index")+
  scale_x_discrete(labels=c("male" = "M", "female" = "F"))+
  theme(legend.position = "none", aspect.ratio = 3.5)
p3

library(gridExtra)
grid.arrange(p2, p3, nrow=1)












