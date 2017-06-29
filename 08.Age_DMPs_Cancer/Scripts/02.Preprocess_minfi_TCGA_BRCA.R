################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
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
idat <- "/ihome/kjohnson/TCGA_BRCA_BMIQ/IDAT"

# Be sure to have a .csv file of covariates in the same location of the idat files
targets <- read.450k.sheet(idat)

# Reads the 'targets' like a dataframe  
RGsetEx <- read.450k.exp(targets = targets)

# Access the phenotype data associated with an experiment  
pheno <- pData(RGsetEx)  

# Produce a density plot of beta values   
png("TCGA_raw_betas.png")
densityPlot(RGsetEx, main = "Beta", xlab = "Beta")
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
load("/ihome/kjohnson/TCGA_BRCA_BMIQ/annotation-150214.RData")

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
save(Betas_BMIQ, file="/ihome/kjohnson/TCGA_BRCA_BMIQ/TCGA_BRCA_Betas_BMIQ.RData")

# Post FunNorm and BMIQ adjustment
png("TCGA_raw_betas.png")
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
annotation <- read.table("/ihome/kjohnson/TCGA_BRCA_BMIQ/annotationfile.tsv", stringsAsFactors = F, row.names = 2, header = T, sep = "\t", nrows = 500000, comment.char = "")[,-1]
# Retrieve all non-Chen and non-sex specific probes and SNP probes to subset beta file
use <- rownames(annotation[annotation$excludeChen == 0 & annotation$SNPinProbe == 0 & annotation$Sex == 0,])
cat("Remove", (length(cgUse) - length(use)), "due to Chen Probes, SNP Probes, and Sex Probes\n")
# Consider the intersection of the cgs given in the beta file (detection P value filtered) and the excluded cgs
filtered.p  <- intersect(cgUse, use)

# Use the intersect function to eliminate bad probes from FunNorm data.
# This is done using the p-value matrix from the previous steps.
intersection <- intersect(rownames(Betas_BMIQ), filtered.p)
# How many probes remain?  
length(intersection) #390,292

# Subset the FunNorm results
Betas_BMIQ_sub = Betas_BMIQ[intersection,]

#################################################                
# STEP 3: Extract Beta-values              
################################################# 
# Get Beta-Values
Betas_processed_BMIQ <- Betas_BMIQ_sub
# Write .csv file of komen betas. Final file to be used in RefFreeEWAS analyses.
write.csv(Betas_processed_BMIQ, "/ihome/kjohnson/TCGA_BRCA_BMIQ/BRCA_BMIQ_Betas_27July2016.csv")


