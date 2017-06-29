###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Run the RefFreeEWAS package (2.0) on the normal breast tissue (Komen)
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Project working directory
setwd("/ihome/kjohnson/TCGA_Norm_RefFreeEWAS")

# Read in the Komen file (larger file size)
library(readr)
library(RefFreeEWAS)

# Load Betas that have been processed by minfi
TCGA_Betas <- read_csv("TCGA_BMIQ_Betas_25July2016.csv", col_names=TRUE)
rownames(TCGA_Betas) <- TCGA_Betas[,1]
TCGA_Betas <- TCGA_Betas[,-1]
Betas = TCGA_Betas[ ,order(colnames(TCGA_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("TCGA_Normal_Breast_manifest_26Feb2016.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)

# Keep all samples in analysis 
covars_TCGA <- covars
Betas_TCGA <- data.matrix(Betas)
# Discount double check
all(names(Betas_TCGA)== covars_TCGA$Match_IDs)

# Choose only the top 10,000 rows (i.e. those with the largest var)
TCGA_Var <- apply(Betas_TCGA, 1, var)
Betas_TCGA_Var = Betas_TCGA[order(TCGA_Var, decreasing=T), ]
Y_shortened = Betas_TCGA_Var[1:10000, ]

# Step 1 - 2: Alternate fixing Mu and Omega by iterating from 2 to 10 (Kmax) cell types
TCGA_RefFree_Array <- RefFreeCellMixArray(Y_shortened, Klist=2:10, iters=25)

# Step 3: Bootstrap method for determining the optimal number of Classes K
RefFRee_TCGA_Boots = RefFreeCellMixArrayDevianceBoots(TCGA_RefFree_Array, Y_shortened, R=1000, bootstrapIterations=10)

# Step 4: Identify the minimum value and store it to be run in short on the local
Compare_K = apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.25)
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.25))

# Save Results
save(list=c("TCGA_RefFree_Array", "RefFRee_TCGA_Boots", "Compare_K"), 
     file="RefFree-TCGA-26July2016.RData", compress=TRUE)