###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Run the RefFreeEWAS package (2.0) on the normal breast tissue (NDRI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Package
library(RefFreeEWAS)
library(data.table)

# Working directory for validation data
setwd("/Users/kevinjohnson/Komen_Normal_5mC/03.NDRI_Validation/")

# Load Betas that have been processed by minfi
NDRI_Betas <- data.frame(fread("Files/NDRI_BMIQ_Betas_25July2016.csv"), row.names=1)
colnames(NDRI_Betas) <- substring(colnames(NDRI_Betas), 2, length(colnames(NDRI_Betas)))
Betas = NDRI_Betas[ ,order(colnames(NDRI_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("Files/metaData-breast-160209-BSonly.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]
covars$Sample_Age <- as.numeric(covars$Sample_Age)

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)
NDRI_Betas <- as.matrix(Betas)

# Choose only the top 10,000 rows (i.e. those with the largest var)
NDRI_Var <- apply(NDRI_Betas, 1, var)
Betas_NDRI_Var = NDRI_Betas[order(NDRI_Var, decreasing=T), ]
Y_shortened = Betas_NDRI_Var[1:10000, ]

# Step 1 - 2: Alternate fixing Mu and Omega by iterating from 2 to 10 (Kmax) cell types
NDRI_RefFree_Array <- RefFreeCellMixArray(Y_shortened, Klist=2:10, iters=25)

# Step 3: Bootstrap method for determining the optimal number of Classes K
RefFRee_NDRI_Boots = RefFreeCellMixArrayDevianceBoots(NDRI_RefFree_Array, Y_shortened, R=1000, bootstrapIterations=10)

# Step 4: Identify the minimum value and store it to be run in short on the local
Compare_K = apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.25)
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.25))

# Save Results
save(list=c("NDRI_RefFree_Array", "RefFRee_NDRI_Boots", "Compare_K"), 
     file="Files/RefFree2-NDRI-25July2016.RData", compress=TRUE)


