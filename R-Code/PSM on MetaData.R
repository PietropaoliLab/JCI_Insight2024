library(MatchIt)
library(tableone)
library(phyloseq)
library(here)
load("AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PlyloSeq Object for analysys.RData")
MetaData <- data.frame(ps@sam_data)





perio <- subset(MetaData, Periodontitis == "Periodontitis")
healthy <- subset(MetaData, Periodontitis == "Healthy")
  
## List numerically coded categorical variables
factorVars <- c("Sex","Site")

## Create a variable list. Use dput(names(pbc))
vars <- c("Sex", "AgeYrs", "Periodontitis", "Smoking", "Race", "Site", "Geo", "BioProject")


## Create Table 1 stratified by trt (omit strata argument for overall table)
tableOne <- CreateTableOne(vars = vars, strata = "Sex", data = MetaData, factorVars = factorVars)

print(tableOne,  showAllLevels = TRUE, noSpaces = TRUE, quote = TRUE)



# In questo PSM devo invertire 0 ed 1 senno non funziona
perio$Male <- ifelse(perio$Sex == "male", 1, 0)
healthy$Male <- ifelse(healthy$Sex == "male", 1, 0)


# Choose subjects according to propensity score
set.seed(101)
match_model <- matchit(Male ~ AgeYrs,
                       data = perio,
                       method = "nearest",
                       ratio = 1,
                       distance = "GLM",
                       link = "logit")
match_data <- match.data(match_model)
colnames(match_data)[ncol(match_data)] <- "propensity_weights"
psm_perio <- perio[perio$Run %in% match_data$Run,]


summary(match_model)
plot(match_model)


## Create Table 1 stratified by trt (omit strata argument for overall table)
tableOne <- CreateTableOne(vars = vars, strata = "Sex", data = psm_perio, factorVars = factorVars)
print(tableOne,  showAllLevels = TRUE, noSpaces = TRUE, quote = TRUE)

MyRun <- psm_perio$Run

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PlyloSeq Object for analysys.RData")
ps
psm_perio <- subset_samples(ps, Run %in% MyRun)
psm_perio
save(psm_perio, file = "Microbiome/Gender features of periodontal microbiome/Data/PSM Periodontitis.RData")






# Choose subjects according to propensity score
set.seed(101)
match_model <- matchit(Male ~ AgeYrs,
                       data = healthy,
                       method = "nearest",
                       ratio = 1,
                       caliper = .25,
                       distance = "GLM",
                       link = "logit")
match_data <- match.data(match_model)
colnames(match_data)[ncol(match_data)] <- "propensity_weights"
psm_healthy <- healthy[healthy$Run %in% match_data$Run,]


## Create Table 1 stratified by trt (omit strata argument for overall table)
tableOne <- CreateTableOne(vars = vars, strata = "Sex", data = psm_healthy, factorVars = factorVars)
print(tableOne,  showAllLevels = TRUE, noSpaces = TRUE, quote = TRUE)

MyRun <- psm_healthy$Run

load("~/Google Drive/01- UNIVERSITA/AQ Projects/Microbiome/Gender features of periodontal microbiome/Data/PlyloSeq Object for analysys.RData")
psm_healthy <- subset_samples(ps, Run %in% MyRun)
ps
psm_healthy

save(psm_healthy, file = "Microbiome/Gender features of periodontal microbiome/Data/PSM Healthy.RData")



