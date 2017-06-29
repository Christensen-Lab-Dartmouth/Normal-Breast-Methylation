###########################
# Breast cancer risk factors are associated with DNA methylation in non-diseased breast tissue independent of cell type
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Run the RefFreeEWAS2.0 from Houseman et al 2016 on normal breast tissue (Komen). 
# Computationally intensive so it is advised to run this script in an HPC environment
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Necessary packages
library(data.table)
library(RefFreeEWAS)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load Betas that have been processed by minfi
Komen_Betas <- data.frame(fread("01.Preprocessing/Data/Komen_BMIQ_Betas_17May2016.csv"), row.names=1)
colnames(Komen_Betas) <- substring(colnames(Komen_Betas), 2, length(colnames(Komen_Betas)))
Betas = Komen_Betas[ ,order(colnames(Komen_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("01.Preprocessing/Data/Komen_Manifest_20July2015.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)

# Keep all samples in analysis 
covars_Komen <- covars
Betas_Komen <- data.matrix(Betas) # data needs to be formatted as matrix
# Discount double check of matching IDs
all(names(Betas_Komen)==covars_Komen$Match_IDs)

# Choose only the top 10,000 CpGs (i.e. those with the largest variance). These are likely to discriminate cell-types.
Komen_Var <- apply(Betas_Komen, 1, var)
Betas_Komen_Var = Betas_Komen[order(Komen_Var, decreasing=T), ]
Y_shortened = Betas_Komen_Var[1:10000, ] # Y represents DNA methylation matrix in model: Y = Mu*Omega

# Step 0: You can start with an intial estimate of M (Optional)
# Komen_Initialized = RefFreeCellMixInitialize(Betas_Komen, K=6)

# Step 1 - 2: Alternate fixing Mu and Omega by iterating from 2 to 10 (Kmax) cell types
Komen_RefFree_Array <- RefFreeCellMixArray(Y_shortened, Klist=2:10, iters=25)
 
# Step 3: Bootstrap method for determining the optimal number of Classes K
RefFree_Komen_Boots = RefFreeCellMixArrayDevianceBoots(Komen_RefFree_Array, Y_shortened, R=1000, bootstrapIterations=10)

# Save Results
save(list=c("Komen_RefFree_Array", "RefFRee_Komen_Boots"), 
     file="RefFree2-Komen-18May2016.RData", compress=TRUE)


