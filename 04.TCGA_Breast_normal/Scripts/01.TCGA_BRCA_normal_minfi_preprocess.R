################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Pre-process the adjacent-to-tumor normal breast TCGA 450K data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~
################################
#Load Libraries
################################
library(minfi) # 1.16.1
library(IlluminaHumanMethylation450kmanifest) 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(wateRmelon)


##########################################
#Step 1: Load IDAT files of interest
##########################################
# Location of IDAT files for study (Place the IDAT files in a designated file location)
idat <- "/Users/kevinjohnson/Komen_Normal_5mC/04.TCGA_Breast_Normal/TCGA_IDATS/IDATS"

#Be sure to have a .csv file of covariates in the same location of the idat files
targets <- read.450k.sheet(idat)

#Reads the 'targets' like a dataframe  
RGsetEx <- read.450k.exp(targets = targets)

#Access the phenotype data associated with an experiment  
pheno <- pData(RGsetEx)  

#Produce a density plot of beta values       
densityPlot(RGsetEx, sampGroups= pheno$ajcc_pathologic_tumor_stage, main = "Beta", xlab = "Beta")
mdsPlot(RGsetEx)

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

TCGA_Betas_BMIQ <- Betas

# Application of BMIQ is computationally slow
save(TCGA_Betas_BMIQ, file="/Users/kevinjohnson/Komen_Normal_5mC/04.TCGA_Breast_Normal/Files/TCGA_Betas_BMIQ.RData")

# Post FunNorm and BMIQ adjustment
densityPlot(TCGA_Betas_BMIQ, main = "TCGA_BMIQ", xlab = "Beta")

#Filtering Probes:
#Probes with a detection p-value above signifigance (0.00001) should approached with caution
detect.p = detectionP(RGsetEx, type = "m+u")
#Begin Filtering Steps at a Detection P value of 0.00001
failed <- detect.p > 0.00001
#The detection p value must be above 0.00001 for 25% of all samples to remove the probe
cg.f <- rowMeans(failed) > 0.25
#Subset the cgs used 
cgUse <- rownames(detect.p)[!cg.f] #Remove 3,112 cgs
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
intersection <- intersect(rownames(TCGA_Betas_BMIQ), filtered.p)
# How many probes remain?  
length(intersection) # 388,596

# Subset the FunNorm results
TCGA_Betas_BMIQ_sub = TCGA_Betas_BMIQ[intersection,]
#################################################                
# STEP 3: Extract Beta-values              
################################################# 
# Get Beta-Values
TCGA_Betas_processed_BMIQ <- TCGA_Betas_BMIQ_sub

# Write .csv file of NDRI betas. Final file to be used in RefFreeEWAS analyses.
write.csv(TCGA_Betas_processed_BMIQ, "/Users/kevinjohnson/Komen_Normal_5mC/04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Betas_25July2016.csv")

