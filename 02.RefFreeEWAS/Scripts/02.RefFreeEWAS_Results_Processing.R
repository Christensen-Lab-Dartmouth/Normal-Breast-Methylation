###########################
# Breast cancer risk factors are associated with DNA methylation in non-diseased breast tissue independent of cell type
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Process the results from the RefFreeEWAS package (2.0) that have presumably been run on a computing cluster
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)
library(RefFreeEWAS)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load in results for 1,000 bootstraps generated in "01.NMF_RefFreeEWAS.R"
load("02.RefFreeEWAS/Files/RefFree-2-18May2016.RData")

# Step 3: Bootstrap method for determining the optimal number of Classes K
RefFRee_Komen_Boots

# Step 4: Identify the minimum value and store it to be run in short on the local
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, median)) # top number is the putative cell type number, bottom number is element of vector indexed
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, mean))
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, mean, trim=0.05))
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, mean, trim=0.15))
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, mean, trim=0.25))
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, mean, trim=0.35))
which.min(apply(RefFRee_Komen_Boots[-1, ], 2, mean, trim=0.45)) # Here, all minimized values were 6 cell-types

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
Betas_Komen <- data.matrix(Betas)

# Sanity check for arrangement of samples
all(names(Betas_Komen)==covars_Komen$Match_IDs)

# Y Matrix of m CpGs x n Subjects of DNA methylation Beta-values
Keq_BMIQ <- RefFreeCellMix(Betas_Komen, mu0=RefFreeCellMixInitialize(Betas_Komen, K=6))

# Examine decomposed matrices of DNAme (i.e., Mu) and cellular proportions (i.e., Omega)
head(Keq_BMIQ$Mu)
head(Keq_BMIQ$Omega)

# Visualize the distribution of cellular proportions
boxplot(Keq_BMIQ$Omega)

# Set any negative cellular proportions to ZERO
Keq_BMIQ$Omega[Keq_BMIQ$Omega < 0] = 0

# All sample IDs match?
all(rownames(Keq_BMIQ$Omega)==covars_Komen$Match_IDs)
Komen_Cell_Proportions = cbind(covars_Komen, Keq_BMIQ$Omega)
write.csv(Komen_Cell_Proportions, file = "/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Files/Komen_Cell_Proportions_BMIQ.csv")

# Use custom heatmap.3 function
source("/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Scripts/07.make_heatmaps.R")

# Age at donation  
colorAge <- hsv(0.01, covars_Komen$Age/max(covars_Komen$Age), 1)

# BMI at donation 
colorBMI <- hsv(0.1, covars_Komen$BMI/max(covars_Komen$BMI), 1)
rowcolors <- as.matrix(t(rbind(colorAge, colorBMI)))
colnames(rowcolors) <- c("Subject Age", "Subject BMI")

# Define the color palette
heatcols <- colorRampPalette(c("white","purple"))(n = 1000)
ppi <- 300
png("02.RefFreeEWAS/Figures/Cell_Proportions_BMIQ_Komen.png", width=7*ppi, height=7*ppi, res=ppi)
heatmap.3(Keq_BMIQ$Omega, col=heatcols, KeyValueName="Cell Proportion", keysize=1.25, dendrogram='col', main = "", 
          margins=c(5,7), xlab="K (Putative Cell-type)", labRow=FALSE)
dev.off()


