###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Permutation testing for associations between Metadata and Omega in NDRI breast samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Packages
library(RefFreeEWAS)
library(data.table)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC")

# Load results from RefFree (K=2:6)
load("03.NDRI_Validation/Files/RefFree2-NDRI-25July2016.RData")
NDRI_RefFree_Array
lapply(NDRI_RefFree_Array, summary)

# Step 3: Bootstrap method for determining the optimal number of Classes K
RefFRee_NDRI_Boots

# Step 4: Identify the minimum value and store it to be run in short on the local
apply(RefFRee_NDRI_Boots[-1, ], 2, mean)
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, median))
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean))
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.05))
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.15))
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.25))
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.35))
which.min(apply(RefFRee_NDRI_Boots[-1, ], 2, mean, trim=0.45))

# 2 putative cell-types appears to be optimal in the NDRI
# Load Betas that have been processed by minfi (No BMIQ)
NDRI_Betas <- data.frame(fread("03.NDRI_Validation/Files/NDRI_BMIQ_Betas_25July2016.csv"), row.names=1)
colnames(NDRI_Betas) <- substring(colnames(NDRI_Betas), 2, length(colnames(NDRI_Betas)))
Betas = NDRI_Betas[ ,order(colnames(NDRI_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("03.NDRI_Validation/Files/metaData-breast-160209-BSonly.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]
covars$Sample_Age <- as.numeric(covars$Sample_Age)

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)
NDRI_Betas <- as.matrix(Betas)

# Discount double check
all(names(NDRI_Betas)== covars$Match_IDs)

# Initialize RefFreeCellMix with the identified optimal cell type number
NDRI_Keq <- RefFreeCellMix(NDRI_Betas, mu0=RefFreeCellMixInitialize(NDRI_Betas, K=2))

head(NDRI_Keq$Omega)
NDRI_Keq$Omega[NDRI_Keq$Omega < 0] = 0
NDRI_Omega <- NDRI_Keq$Omega

# Ensure that IDs match one another
all(rownames(NDRI_Omega)==covars$Match_IDs)

# Exportable cellular proportions (K hat)
NDRI_Cell_Proportions = cbind(covars, NDRI_Omega)
write.csv(NDRI_Cell_Proportions, file="03.NDRI_Validation/Files/NDRI_BMIQ_Cell_Proportions.csv")

# Generate heatmap for cellular proportions
# Use heatmap.3 function
source("/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Scripts/07.make_heatmaps.R")

#Age at death (medical record)
colorAge <- hsv(0.01, covars$Sample_Age/max(covars$Sample_Age), 1)

#BMI at death (medical record)
colorBMI <- hsv(0.1, covars$Sample_BMI/max(covars$Sample_BMI), 1)

rowcolors <- as.matrix(t(rbind(colorAge, colorBMI)))
colnames(rowcolors) <- c("Subject Age", "Subject BMI")

#Define the color palette
heatcols <- colorRampPalette(c("white","purple"))(n = 1000)
ppi = 300
png("03.NDRI_Validation/Figures/Cell_Proportions_BMIQ_NDRI.png", width=7*ppi, height=7*ppi, res=ppi)
heatmap.3(NDRI_Omega, col=heatcols, KeyValueName="Cell Proportion", keysize=1.25, dendrogram='col', main = "", 
          margins=c(5,7), xlab="K (Putative Cell-type)", labRow=FALSE)
dev.off()
