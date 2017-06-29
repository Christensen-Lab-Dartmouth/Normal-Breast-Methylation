###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify optimal Omega in TCGA adjacent normal breast samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Packages
library(RefFreeEWAS)
library(data.table)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load results from RefFree (K=2:10)
load("04.TCGA_Breast_Normal/Files/RefFree-TCGA-26July2016.RData")
TCGA_RefFree_Array
lapply(TCGA_RefFree_Array, summary)

# Step 3: Bootstrap method for determining the optimal number of Classes K
RefFRee_TCGA_Boots

# Step 4: Identify the minimum value and store it to be run in short on the local
apply(RefFRee_TCGA_Boots[-1, ], 2, mean)
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, median))
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean))
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.05))
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.15))
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.25))
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.35))
which.min(apply(RefFRee_TCGA_Boots[-1, ], 2, mean, trim=0.45))

# 7 putative cell-types appears to be optimal in the TCGA normal-adjacent breast

# Load Betas that have been processed by minfi
TCGA_Betas <- data.frame(fread("04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Betas_25July2016.csv"), row.names=1)
colnames(TCGA_Betas) <- substring(colnames(TCGA_Betas), 2, length(colnames(TCGA_Betas)))
Betas = TCGA_Betas[ ,order(colnames(TCGA_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("04.TCGA_Breast_Normal/Files/2016-08-04_TCGA_Normal_Breast_manifest.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)

# Keep all samples in analysis 
covars_TCGA <- covars
Betas_TCGA <- data.matrix(Betas)
# Discount double check
all(names(Betas_TCGA)== covars_TCGA$Match_IDs)

# Initialize RefFreeCellMix with the identified optimal cell type number
TCGA_Norm_Keq <- RefFreeCellMix(Betas_TCGA, mu0=RefFreeCellMixInitialize(Betas_TCGA, K = 7))

head(TCGA_Norm_Keq$Omega)
TCGA_Norm_Keq$Omega[TCGA_Norm_Keq$Omega < 0] = 0
TCGA_Omega <- TCGA_Norm_Keq$Omega

# Ensure that IDs match one another
all(rownames(TCGA_Omega)==covars$Match_IDs)

# Exportable cellular proportions (K hat)
TCGA_Cell_Proportions = cbind(covars, TCGA_Omega)
write.csv(TCGA_Cell_Proportions, file="04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Cell_Proportions.csv")
TCGA_Cell_Proportions <- read.csv("04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Cell_Proportions.csv", header=T)
TCGA_Omega <- as.matrix(TCGA_Cell_Proportions[ ,11:17])
colnames(TCGA_Omega) <- as.character(seq(1,7))
  
# Generate heatmap for cellular proportions
# Use heatmap.3 function
source("/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Scripts/07.make_heatmaps.R")

#Age at death (medical record)
colorAge <- hsv(0.01, covars$Age/max(covars$Age), 1)

rowcolors <- as.matrix(t(colorAge))
names(rowcolors) <- c("Age")

#Define the color palette
heatcols <- colorRampPalette(c("white","purple"))(n = 1000)
ppi = 300
png("04.TCGA_Breast_Normal/Figures/TCGA_BRCA[n]_Cellular_Proportions.png", width=7*ppi, height=7*ppi, res=ppi)
heatmap.3(TCGA_Omega, labRow= "", col=heatcols, KeyValueName="Cell Proportion", dendrogram='column', trace = "none",
          main = "TCGA Normal", margins=c(5,7), xlab="K (Putative Cell-Type)")
dev.off()
