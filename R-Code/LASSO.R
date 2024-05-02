library(tidyverse)
library(SIAMCAT)
library(ggplot2)
library(phyloseq)
library(microbiome)
# MacOs
 load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")

# windows 
# load('G:/My Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData')


# https://www.bioconductor.org/packages/devel/bioc/vignettes/SIAMCAT/inst/doc/SIAMCAT_meta.html

sc.list <- list()
ps <- psm_perio

sites <- unique(ps@sam_data$Site)

for (s in sites) {
  s="Saliva"
  ps <- subset_samples(ps, Site == s)
  
sample_data(ps) <- sample_data(ps)[,c('Run','BioProject', 'Sex')]

# Convert in relative abundance
ps <- transform(ps, transform = "compositional")

meta.all <- data.frame(sample_data(ps))
feat <- otu_table(ps)
feat <- t(as.matrix(feat))


meta.ind <- meta.all %>% 
              group_by(Run) %>% 
              ungroup()

# create tibble to store all the predictions
auroc.all <- tibble(study.train=character(0), 
                    study.test=character(0),
                    AUC=double(0))


# and a list to save the trained SIAMCAT objects

  datasets <- unique(meta.all$BioProject)
  

for (i in datasets){
  i="PRJNA321534" 
  # restrict to a single study
  meta.train <- meta.all %>% 
    filter(BioProject==i) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Run
  
  
  # run the constructor function
  sc.obj.train <- siamcat(feat=feat, meta=meta.train, 
                          label='Sex', case='female')


  # normalize features
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
                                     norm.param=list(log.n0=1e-05, sd.min.q=0),feature.type = 'original')
  # Create data split
  sc.obj.train <- create.data.split(sc.obj.train,
                                    num.folds = 10, num.resample = 3)
  # train LASSO model
  sc.obj.train <- train.model(sc.obj.train, method='lasso')
  
  ## apply trained models to other datasets
  
  # loop through datasets again
  for (i2 in datasets){
    i2=i 
    if (i == i2){
      # make and evaluate cross-validation predictions (same dataset)
      sc.obj.train <- make.predictions(sc.obj.train)
      sc.obj.train <- evaluate.predictions(sc.obj.train)
      auroc.all <- auroc.all %>% 
        add_row(study.train=i, study.test=i, 
                AUC=eval_data(sc.obj.train)$auroc %>% as.double())
    } else {
      # make and evaluate on the external datasets
      # use meta.ind here, since we want only one sample per subject!
      meta.test <- meta.ind %>% 
        filter(BioProject==i2) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Run
      
      sc.obj.test <- siamcat(feat=feat, meta=meta.test,
                             label='Sex', case='female')
      # make holdout predictions
      sc.obj.test <- make.predictions(sc.obj.train, 
                                      siamcat.holdout = sc.obj.test)
      sc.obj.test <- evaluate.predictions(sc.obj.test)
      auroc.all <- auroc.all %>% 
        add_row(study.train=i, study.test=i2, Site = s,
                AUC=eval_data(sc.obj.test)$auroc %>% as.double())
    }
  }
  # save the trained model
  sc.list[[paste0(i,"_", s)]] <- sc.obj.train
}

}


test.average <- auroc.all %>% 
  filter(study.train!=study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC=mean(AUC), .groups='drop') %>% 
  mutate(study.train="Average")




# combine AUROC values with test average
df <- bind_rows(auroc.all, test.average) %>% 
  # highlight cross validation versus transfer results
  mutate(CV=study.train == study.test) %>%
  # for facetting later
  mutate(split=case_when(study.train=='Average'~'Average', TRUE~'none')) %>% 
  mutate(split=factor(split, levels = c('none', 'Average'))) %>% 
  # convert to factor to enforce ordering
  mutate(study.train=factor(study.train, levels=c(datasets, 'Average'))) %>% 
  mutate(study.test=factor(study.test, 
                           levels=c(rev(datasets),'Average'))) 

p1 <- ggplot(df, aes(y=study.test, x=study.train, fill=AUC, size=CV, color=CV)) +
  geom_tile() + theme_minimal() +
  # text in tiles
  geom_text(aes_string(label="format(AUC, digits=2)"), 
            col='white', size=2)+
  # color scheme
  scale_fill_distiller(type = 'seq', palette = 'BuGn', limits=c(0.5, 1)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Training Set') + ylab('Test Set') + 
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 1), guide=FALSE) + 
  facet_grid(~split, scales = 'free', space = 'free')

p1











