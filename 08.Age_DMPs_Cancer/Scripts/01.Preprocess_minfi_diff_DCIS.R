################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
#################################
# ~~~~~~~~~~~~~~~~~~
# The script will curate useable samples and CpG sites and preprocess with Funnorm
# Using the package minfi as introducted by Aryee et al 2014 Workflow as described in the manuscript
# ~~~~~~~~~~~~~~~~~~
# To install minfi first install bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite()
################################
# Load Libraries
################################
library(minfi) # minfi_1.16.1
library(IlluminaHumanMethylation450kmanifest) 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(wateRmelon) 

##########################################
# Step 1: Load IDAT files of interest
##########################################
# Location of IDAT files for study (Place the IDAT files in a designated file location)
idat <- "/Users/kevinjohnson/Dropbox (Christensen Lab)/DCIS_Study_2013/Dec_DCIS_minfi"

# Be sure to have a .csv file of covariates in the same location of the idat files
targets <- read.450k.sheet(idat)

# Reads the 'targets' like a dataframe  
RGsetEx <- read.450k.exp(targets = targets)

# Access the phenotype data associated with an experiment  
pheno <- pData(RGsetEx)  

# Produce a density plot of beta values   
png("/Users/kevinjohnson/Komen_Normal_5mC/07.Age_DMPs_Cancer/Figures/NHMN_raw_betas.png")
densityPlot(RGsetEx, sampGroups = pheno$Sample_Group, main = "Beta", xlab = "Beta")
dev.off()

#################################################                
# STEP 2: Convert intensities to methylation signal and normalize              
################################################# 
# Functional Normalization is a between-array normalization. Removes unwanted variation by regressing
# out variability explained by the control probes present on the array. 
Mset_FunNorm <- preprocessFunnorm(RGsetEx, nPCs=2, sex = pheno$gender)

# Extract beta-values from the Funnorm methylation set
Betas <- getBeta(Mset_FunNorm)

# The BMIQ function is located in the wateRmelon package
# the long annotation file for the 450K array was prepared by EA Houseman to identify probe-type
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Files/annotation-150214.RData")

# Identify probe types
infType <- longAnnot$Type
names(infType) <- longAnnot$TargetID
infType <- infType[rownames(Betas)]

# Change naming scheme 
ttp <- ifelse(infType=="I",1,2)

n <- dim(Betas)[2]
for(i in 1:n){
  cat(i,"\n")
  try({
    flag <- which(!is.na(Betas[,i]))
    bq <- BMIQ(Betas[flag,i],ttp[flag])
    Betas[flag,i] <- bq$nbeta
  })
}

Betas_BMIQ <- Betas

# Application of BMIQ is computationally slow
save(Betas_BMIQ, file="/Users/kevinjohnson/Komen_Normal_5mC/07.Age_DMPs_Cancer/Files/NHMN_DCIS_Betas_BMIQ.RData")

# Post FunNorm and BMIQ adjustment
png("/Users/kevinjohnson/Komen_Normal_5mC/07.Age_DMPs_Cancer/Figures/NHMN_DCIS_Betas_BMIQ.png")
densityPlot(Betas_BMIQ, main = "Beta_BMIQ", xlab = "Beta")
dev.off()

# Filtering Probes:
# Probes with a detection p-value above signifigance (0.00001) should approached with caution
detect.p = detectionP(RGsetEx, type = "m+u")
# Begin Filtering Steps at a Detection P value of 0.00001
failed <- detect.p > 0.00001
# The detection p value must be above 0.00001 for 25% of all samples to remove the probe
cg.f <- rowMeans(failed) > 0.25
# Subset the cgs used 
cgUse <- rownames(detect.p)[!cg.f] #Remove 489 cgs
cat("Remove", (nrow(detect.p) - length(cgUse)), "due to high detection p values\n")

# Remove Chen Probes and Sex Specific Probes
# Load Annotation File
annotation <- read.table("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Files/annotationfile.tsv", stringsAsFactors = F, row.names = 2, header = T, sep = "\t", nrows = 500000, comment.char = "")[,-1]
# Retrieve all non-Chen and non-sex specific probes and SNP probes to subset beta file
use <- rownames(annotation[annotation$excludeChen == 0 & annotation$SNPinProbe == 0 & annotation$Sex == 0,])
cat("Remove", (length(cgUse) - length(use)), "due to Chen Probes, SNP Probes, and Sex Probes\n")
# Consider the intersection of the cgs given in the beta file (detection P value filtered) and the excluded cgs
filtered.p  <- intersect(cgUse, use)

# Use the intersect function to eliminate bad probes from FunNorm data.
# This is done using the p-value matrix from the previous steps.
intersection <- intersect(rownames(Betas_BMIQ), filtered.p)
# How many probes remain?  
length(intersection) # 371762

# Subset the FunNorm results
Betas_BMIQ_sub = Betas_BMIQ[intersection,]

#################################################                
# STEP 3: Extract Beta-values              
################################################# 
# Get Beta-Values
Betas_processed_BMIQ <- Betas_BMIQ_sub
# Write .csv file of komen betas. Final file to be used in RefFreeEWAS analyses.
write.csv(Betas_processed_BMIQ, "/Users/kevinjohnson/Komen_Normal_5mC/07.Age_DMPs_Cancer/Files/DCIS_BMIQ_Betas.csv")

#################################################                
# STEP 4: apply limma fit               
################################################# 
library(limma)
library(qvalue)
library(readr)
Betas_processed_BMIQ <- read_csv("/Users/kevinjohnson/Komen_Normal_5mC/07.Age_DMPs_Cancer/Files/DCIS_BMIQ_Betas.csv", col_names=T)
rownames(Betas_processed_BMIQ) <- Betas_processed_BMIQ[,1]
Betas_processed_BMIQ <- Betas_processed_BMIQ[,-1]

# Examine all CpGs
DCIS_Betas <- as.matrix(Betas_processed_BMIQ[, ]) 
# Discount double check
all(names(DCIS_Betas)==rownames(pheno))
XX <- model.matrix(~Sample_Group + agepr, data= pheno)

# Convert Komen beta-values to M-values
Betas_DCISM <- ifelse(DCIS_Betas>=1,1-1E-6,ifelse(DCIS_Betas<=0,1E-6, DCIS_Betas))
Betas_DCISM <- log(Betas_DCISM)-log(1-Betas_DCISM)

all(colnames(Betas_DCISM)==rownames(XX))
lf_DCIS <-  eBayes(lmFit(Betas_DCISM, XX))

setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
overlap_CpGs <- Age_DMPs[which(Age_DMPs%in%rownames(lf_DCIS$p.value))]

hist(lf_DCIS$p.value[ , 'Sample_GroupDCIS'])
hist(lf_DCIS$p.value[overlap_CpGs , 'Sample_GroupDCIS'])

# Tally the number of CpGs that are below the q-value 0.01 threshold
q_dcis <- qvalue(lf_DCIS$p.value[ , 'Sample_GroupDCIS'], fdr.level=0.01)
table(q_dcis$significant) # 31,186


# Sample for significant difference (P-values). Not a statistically appropriate approach, just exploration
set.seed(123)
random_sites <- sample(lf_DCIS$p.value[ , 'Sample_GroupDCIS'], 775)
age_related_sites <- lf_DCIS$p.value[overlap_CpGs , 'Sample_GroupDCIS']
wilcox.test(age_related_sites,random_sites) # 9.526e-05

# Sample for significant difference (Coefficients). Assessment of how large the effects are in age-related vs. random
set.seed(123)
random_sites <- abs(sample(lf_DCIS$coefficients[ , 'Sample_GroupDCIS'], 775))
age_related_sites <- abs(lf_DCIS$coefficients[overlap_CpGs , 'Sample_GroupDCIS'])
wilcox.test(age_related_sites, random_sites) # 2.134e-09

# Extract the median P-value for comparing the absolute values of coefficients from age-related
# and randomly selected CpGs from the DCIS normal-cancer comparisons
wilcotestgenerator <- function(n) {
  # note that here we have a "age-related" group and a 
  # randomly selected group. The age-related remains the same (n=775)
  age_related_sites <- abs(lf_DCIS$coefficients[overlap_CpGs , 'Sample_GroupDCIS'])
  random_sites <- abs(sample(lf_DCIS$coefficients[ , 'Sample_GroupDCIS'], n))
  wilco_stat <- wilcox.test(age_related_sites, random_sites) 
  return(wilco_stat$p.value)
}
wilco_tests <- replicate(1000, wilcotestgenerator(775))
summary(wilco_tests)

# Almost all of the tests are highly statistically significant
hist(-log10(wilco_tests), xlim=c(0, max(-log10(wilco_tests))))
abline(v= -log10(.05), col="red", lwd=2)





