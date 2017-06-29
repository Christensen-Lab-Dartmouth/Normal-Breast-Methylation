###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# EpiTOC values for the 100 Komen samples. Test for relation with risk factors.
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)
# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load Betas that have been processed by minfi pipeline and then BMIQ (type II probe adjustment)
Komen_Betas <- data.frame(fread("01.Preprocessing/Data/Komen_BMIQ_Betas_17May2016.csv"), row.names=1)
colnames(Komen_Betas) <- substring(colnames(Komen_Betas), 2, length(colnames(Komen_Betas)))
Betas = Komen_Betas[ ,order(colnames(Komen_Betas), decreasing=T)]

# Load covariate file for Komen
covariates <- read.csv("01.Preprocessing/Files/Komen_Manifest_20July2015.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)

# Keep all samples in analysis. Need to match on IDAT file ID
covars_Komen <- covars
rownames(covars_Komen) <- covars_Komen$Match_IDs
Betas_Komen <- data.matrix(Betas)
# Are all the samples in both matrices aligned.
all(names(Betas_Komen)==covars_Komen$Match_IDs)

# EpiTOC CpGs. It is also loaded by the function provided.
load("/Users/kevinjohnson/Komen_Normal_5mC/06.Epigenetic_Clock/Files/EpiTOC_CpGs.rdata")
# Inspect how many age-related breast CpGs overlap with EpiTOC CpGs
load("/Users/kevinjohnson/Komen_Normal_5mC/05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
sum(Age_DMPs%in%epiTOCcpgs.v) # 3 out of 385 MEASURED CpGs

#### DESCRIPTION
#### An R-function to estimate the pcgtAge-score under the EpiTOC model. Only required input argument is a data matrix of BMIQ (or other type-2 probe correction) normalised DNAm data. Output is the pcgtAge-score for each column (sample) in the input data matrix.

#### INPUT:
#### data.m: DNAm data beta-valued matrix with rownames labeling CpGs and columns labeling samples for which a pcgtAge-score is desired
#### mode: this is an optional parameter, which can take on 2 values: "raw" or "z". If the former, function computes the pcgtAge-score as the average methylation over the represented epiTOC loci. If this option is set to "z", then the score is computed as an average of z-scores, reflecting a deviation relative to a specified reference set of samples (as specified by ref.idx)
#### ref.idx: this needs to be specified if mode="z". It is an index vector of column position of samples in data.m which should be treated as a reference.


#### OUPTUT:
#### score: a vector of pcgtAge-scores, one value for each sample of data.m
# Saved rdata file from Yang et al
setwd("/Users/kevinjohnson/Komen_Normal_5mC/06.Epigenetic_Clock/Files/")

EstEpiTOC <- function(data.m,mode="raw",ref.idx=NULL){
  load("EpiTOC_CpGs.rdata");
  common.v <- intersect(rownames(data.m),epiTOCcpgs.v);
  print(paste("Number of represented epiTOC CpGs=",length(common.v),sep=""));
  map.idx <- match(common.v,rownames(data.m));
  if(mode=="raw"){
    score.v <- colMeans(data.m[map.idx,]);
  }
  else if (mode=="z"){
    if(is.null(ref.idx)){
      print("You have selected mode=z, so please specify a reference index vector ref.idx!");
      stop;
      break;
    }
    else {
      sd.v <- apply(data.m[map.idx,ref.idx],1,sd);
      z.m <- (data.m[map.idx,] - rowMeans(data.m[map.idx,ref.idx]))/sd.v;
      score.v <- colMeans(z.m);
    }
    
  }         
  
  return(score.v);
}

# Save output. 356 CpGs out of 385 EpiTOC CpGs were available in Komen data set.
EpiTOCresults_vec <- EstEpiTOC(Betas)
# Provide distribution of EpiTOC values
hist(EpiTOCresults_vec)

# Make sure that EpiTOC didn't realign sample order during analysis
all(names(EpiTOCresults_vec)==covars_Komen$Match_IDs)
# Scatterplot of pcgtAGE vs. Subject Age
png("/Users/kevinjohnson/Komen_Normal_5mC/06.Epigenetic_Clock/Figures/EpiTOC_Komen.png", height=7*300, width=7*300, res=300)
plot(covars_Komen$Age, EpiTOCresults_vec, pch=19, bty='l', ylab="pcgtAge", xlab='Subject Age')
dev.off()

cor.test(covars_Komen$Age, EpiTOCresults_vec) # 0.75

fit = lm(EpiTOCresults_vec~covars_Komen$Pregnant)
summary(fit)

plot(as.factor(covars_Komen$Pregnant), EpiTOCresults_vec)
cor.test(covars_Komen$BMI, EpiTOCresults_vec)


setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
covariates <- read.csv("01.Preprocessing/Files/Komen_Manifest_20July2015.csv", header = T, sep = ",", stringsAsFactors = F,row.names=1)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]
full_covariates <- read.csv("01.Preprocessing/Files/2016_08_08_Full_Covariate_Komen_Samples.csv", header=T, row.names=1)
tmp_covar <- merge(covars, full_covariates, by="row.names")
rownames(tmp_covar) <- tmp_covar$Row.names
tmp_covar$Row.names <- NULL
tmp_covar$Komen_ID <- rownames(tmp_covar)

# Results from submitted 450K array data (unadjusted for BMIQ as Horvath implements normalizatio procedure)
Clock_output <- read.csv("06.Epigenetic_Clock/Files/Betas_Clock_K.output.csv", header=TRUE)

# Properly orient the data
covar_clock <- merge(tmp_covar, Clock_output, by.x="Match_IDs", by.y="X")
rownames(covar_clock) <- covar_clock$Row.names
covar_clock$Row.names <- NULL
horvath_clock <- covar_clock[order(covar_clock$Match_IDs, decreasing = T), ]
all(horvath_clock$Match_IDs==names(EpiTOCresults_vec))

png("/Users/kevinjohnson/Komen_Normal_5mC/06.Epigenetic_Clock/Figures/Horvath_vs_EpiTOC_Komen.png", height=7*300, width=7*300, res=300)
plot(horvath_clock$AgeAccelerationResidual, EpiTOCresults_vec, pch=19, bty='l', xlab="Horvath Age Accel - Residual", ylab="EpiTOC pcgtAge")
dev.off()
cor.test(horvath_clock$AgeAccelerationResidual, EpiTOCresults_vec)

# Change naming scheme for drinks per week
horvath_clock$Drinks.per.Week <- as.factor(horvath_clock$Drinks.per.Week)
horvath_clock$drinking_status <- ifelse(horvath_clock$Drinks.per.Week == " ", 0,  ifelse(horvath_clock$Drinks.per.Week == " < 7", 1, 2) )
Drinks_per_week <- glm(as.numeric(EpiTOCresults_vec) ~ as.factor(horvath_clock$drinking_status),  family= gaussian)
summary(Drinks_per_week) # > 7 drinks per week. P = 0.62
plot(as.factor(horvath_clock$drinking_status), as.numeric(EpiTOCresults_vec))

horvath_clock$Race.x[horvath_clock$Race.x=="AFRNAMER"] = "AFRICAN-AMERICAN"
horvath_clock$Race.x <- relevel(as.factor(horvath_clock$Race.x), ref = "WHITE") # n = 86 white women
Race <- glm(as.numeric(EpiTOCresults_vec)~ as.factor(Race.x), family = gaussian, data = horvath_clock)
summary(Race) # African-American P = 0.021*, Hispanic P = 0.04

horvath_clock$BMI_group <- ifelse(as.numeric(horvath_clock$BMI.x>=25.0), "overweight", "normal") 
table(horvath_clock$BMI_group)
horvath_clock$Drinks.per.Week <- as.factor(horvath_clock$Drinks.per.Week)
horvath_clock$drinking_status <- ifelse(horvath_clock$Drinks.per.Week == " ", 0,  ifelse(horvath_clock$Drinks.per.Week == " < 7", 1, 2) )

Multi_var_full_model <- lm(as.numeric(EpiTOCresults_vec)~as.factor(BMI_group)+as.factor(BLDCAN)+as.factor(Pregnant)+as.factor(Race.x)+as.factor(drinking_status), 
                           data = horvath_clock)
summary(Multi_var_full_model) # B = 4.0, P = 0.035 for African-American AgeAccel. 

# Produce a ggtheme for the background of the plot

Breast_Clock <- ggplot(horvath_clock, aes(x=Age.x, y= EpiTOCresults_vec, color= Race.x , size = 2) )+ geom_point(alpha=1/3) +  ylim(0.025,.10) + xlim(10,100) + 
  geom_abline(intercept= .02, slope=1, colour='black', alpha=1/3)  + theme_bw() 
png("06.Epigenetic_Clock/Figures/EpiTOC-Score-Normal-Figure.png", height=7*300, width=7*300, res=300)
Breast_Clock + annotate('text',x=70, 0.09, label='epiTOC score-Age P = 7.5E-01', alpha=0.8) + xlab("Chronological Age")  + ylab("epiTOC score (DNA Methylation Age)") +
  labs(size="",color="Race") + theme(legend.key = element_blank())
dev.off()




