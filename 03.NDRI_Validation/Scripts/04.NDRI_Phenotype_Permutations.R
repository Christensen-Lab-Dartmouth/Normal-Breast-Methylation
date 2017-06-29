###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Permutation test for determining associations with metadata across classes
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Packages
library(RefFreeEWAS)
library(data.table)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC")

# Load Betas that have been processed by minfi (No BMIQ)
NDRI_Betas <- data.frame(fread("03.NDRI_Validation/Files/NDRI_BMIQ_Betas_25July2016.csv"), row.names=1)
colnames(NDRI_Betas) <- substring(colnames(NDRI_Betas), 2, length(colnames(NDRI_Betas)))
Betas = NDRI_Betas[ ,order(colnames(NDRI_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("03.NDRI_Validation/Files/metaData-breast-160209-BSonly.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T), ]
covars$Sample_Age <- as.numeric(covars$Sample_Age)

# Inspect whether or not the samples are in the same order
all(names(Betas) == covars$Match_IDs)
NDRI_Betas <- as.matrix(Betas)

# Discount double check
all(names(NDRI_Betas) == covars$Match_IDs)

# Build Design Matrix
DESIGNX <- model.matrix(~Sample_Age + Sample_BMI, data = covars)

# Load results from RefFree (K=2:6)
load("03.NDRI_Validation/Files/RefFree2-NDRI-25July2016.RData")
# Examine the structure from the RefFree output
NDRI_RefFree_Array

# Initialize RefFreeCellMix with the identified optimal cell type number
NDRI_Keq <- RefFreeCellMix(NDRI_Betas, mu0=RefFreeCellMixInitialize(NDRI_Betas, K=2))
# save(NDRI_Keq, file="03.NDRI_Validation/Files/NDRI_BMIQ_Keq.RData")
load("03.NDRI_Validation/Files/NDRI_BMIQ_Keq.RData")

#Andy-defined functions for permutation testing
omegArray <- function(omList, xDesign, family = quasibinomial(), 
                      stats=c("Estimate","Pr(>|t|)")){
  nK <- length(omList)
  p <- dim(xDesign)[2]
  nStat <- length(stats)
  results <- list()
  for(i in 1:nK){
    K <- dim(omList[[i]])[2]
    results[[i]] <- array(NA,dim=c(p,nStat,K))
    dimnames(results[[i]]) <- list(colnames(xDesign),stats,1:K) 
    for(k in 1:K){
      omega <- omList[[i]][,k]
      omega <- pmax(pmin(omega,1),0)
      results[[i]][,,k] <- summary(glm(omega~xDesign-1, family=family))$coef[,stats]
    }
  }
  class(results) <- "omegArray"
  results
}

print.omegArray <- function(x){
  cat("OMEGA MODEL SERIES\n")
  cat("Covariates:",paste(dimnames(x[[1]])[[1]],collapse=","),"\n")
  cat("Statistics:",paste(dimnames(x[[1]])[[2]],collapse=","),"\n")
  cat("# Classes:",paste(sapply(x,dim)[3,],collapse=","),"\n")  
}

omegArrayPermutation <- function(omList, xDesign, covariates=2, R=10, ...){
  n <- dim(xDesign)[1]
  est <- omegArray(omList, xDesign, ...)
  perms <- list()
  for(r in 1:R){
    xPermute <- xDesign
    xPermute[,covariates] <- xPermute[sample(1:n,n),covariates]
    perms[[r]] <- lapply(omegArray(omList, xPermute, ...),
                         function(u) {
                           tab <- u[covariates,,,drop=FALSE]
                           dimnames(tab) <- list(colnames(XX)[covariates], dimnames(est[[1]])[[2]], 1:dim(u)[3])
                           tab
                         })
  }
  out <- list(estimate=est, permutations=perms)
  class(out) <- "omegArrayPermutation"
  out
}

print.omegArrayPermutation <- function(x){
  cat("OMEGA MODEL SERIES PERMUTATION\n\nEstimate:\n")
  print(x$estimate)
  cat("\nPermuted covariates: ")
  cat(paste(dimnames(x$permutations[[1]][[1]])[[1]],collapse=","),"\n")
  print(summary(x))
}

minp.functional <- function(coefTableList){
  min(sapply(coefTableList, function(coefTable)min(coefTable[,"Pr(>|t|)",])))
}

make.minp.functional.subclasses <- function(classes){
  function(coefTableList){
    min(sapply(coefTableList[classes], function(coefTable)min(coefTable[,"Pr(>|t|)",])))
  }
}

summary.omegArrayPermutation <- function(x,functional=minp.functional){ 
  permCovs <- dimnames(x$permutations[[1]][[1]])[[1]]
  
  theta0 <- functional(lapply(x$estimate, function(u)u[permCovs,,,drop=FALSE]))
  thetap <- sapply(x$permutations,functional)
  
  out <- mean(thetap<=theta0)
  names(out) <- "p.value"
  out
}

individual.p.values <- function(x){
  nK <- length(x$estimate)
  p <- rep(NA,nK)
  for(k in 1:nK){
    p[k] <- summary(x, functional=make.minp.functional.subclasses(k))
  }
  p
}
#######################################

# NOTE: XBLOCK is a list of indices corresponding to specific (sets of) covariates test
#   where "XX" is the model matrix for the model representing the full model.
### EXAMPLES:
## BREAST TUMORS (#1) (27K)
#XBLOCK <- list(Tumor=2:4,Age=5,Size=6)
#
## BREAST TUMORS (#2) (27K)
#XBLOCK <- list(Tumor=2:6,Age=7)

# Create groups for age by binning into tertiles
summary(covars$Sample_Age)

# If desired, group BMI into two categories split on the median (25.98; near classification for overweight)
summary(covars$Sample_BMI)

## NDRI normal Breast (450K)
XBLOCK <- list(Age = 2, BMI = 3)
XX <- model.matrix(~Sample_Age + Sample_BMI, data = covars)
NDRI_om_List <- sapply(NDRI_RefFree_Array, "[[", 2) 

# Inclusion of High BMI subjects (max BMI = 62, may have to make categorical)
permModels <- lapply(XBLOCK, function(u){
  cat(colnames(XX)[u],"\n")
  omegArrayPermutation(NDRI_om_List[1:5], XX, R=1000, covariates=u)
})
permModels  # Overall P-values for all K
# save(permModels, file="04.NDRI_Validation/Files/Permutation_1000_age_BMI_NDRI.RData")

# Separately by K
sapply(permModels, individual.p.values)

# Over just 1 to 2 (K-hat = 2)
sapply(permModels, summary, functional=make.minp.functional.subclasses(1:2))

# Over just class 1 (should match above "separately by K" for K=1
sapply(permModels, summary, functional=make.minp.functional.subclasses(1))
sapply(permModels, summary, functional=make.minp.functional.subclasses(2))
sapply(permModels, summary, functional=make.minp.functional.subclasses(3))
sapply(permModels, summary, functional=make.minp.functional.subclasses(4))
sapply(permModels, summary, functional=make.minp.functional.subclasses(5))

# Plot for Permutation P-values
png("/Users/kevinjohnson/Komen_Normal_5mC/03.NDRI_Validation/Figures/NDRI_Phenotype_Permutations.png",  width=7*300, height=7*300, res=300)
plot(x=c(2,3,4,5,6), y=-log10(sapply(permModels, individual.p.values)[ , 1]), pch=16, bty="n", col='red',
     ylab='-log10(P-value)', xlab='K (putative cell-types)', ylim=c(0,4.5))
lines(x=c(2,3,4,5,6), y=-log10(sapply(permModels, individual.p.values)[ , 1]), col='red', lwd = 2)
points(x=c(2,3,4,5,6), y=-log10(sapply(permModels, individual.p.values)[ , 2]), pch=16, col='blue')
lines(x=c(2,3,4,5,6), y=-log10(sapply(permModels, individual.p.values)[ , 2]), col='blue', lwd = 2)
text(2, -log10(0.005), "Perm P-val = 0.026", col='red', adj=c(0,-3), cex=0.8)
text(2, -log10(0.048), "Perm P-val = 0.099", col='blue', adj=c(0,-3.5), cex=0.8)
abline(h=-log10(0.05), col = "lightgray", lty = 2, lwd = 2)
text(5.5, -log10(0.04), "P < 0.05 threshold", col="lightgray", cex=0.75)
legend(2, 4, c("Age","BMI"), cex=1, lwd=3, lty = 1, col=c('red', 'blue'), bty='n')
dev.off()

