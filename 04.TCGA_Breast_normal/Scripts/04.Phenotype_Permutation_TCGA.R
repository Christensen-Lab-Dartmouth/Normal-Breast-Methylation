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
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
load("04.TCGA_Breast_Normal/Files/RefFree-TCGA-26July2016.RData")

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
all(names(Betas_TCGA)== covars$Match_IDs)

# Build Design Matrix
DESIGNX <- model.matrix(~Age, data= covars)

# RefFree estimated cell types for TCGA normal breast samples (K-hat = 7)
TCGA_Norm_Keq <- read.csv("04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Cell_Proportions.csv", header=T)

# Andy-defined functions for permutation testing
omegArray <- function(omList, xDesign, family=quasibinomial(), 
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
# EXAMPLES:
# BREAST TUMORS (#1) (27K)
# XBLOCK <- list(Tumor=2:4,Age=5,Size=6)
#
# BREAST TUMORS (#2) (27K)
# XBLOCK <- list(Tumor=2:6,Age=7)
#######################################
# Create groups for age by binning into tertiles
summary(covars$Age)

## NDRI normal Breast (450K)
XBLOCK <- list(Age=2)
XX <- model.matrix(~Age, data = covars)
TCGA_om_List <- sapply(TCGA_RefFree_Array, "[[", 2) 

# Inclusion of all stage subjects
permModels <- lapply(XBLOCK, function(u){
  cat(colnames(XX)[u],"\n")
  omegArrayPermutation(TCGA_om_List[1:9], XX, R=1000, covariates=u)
})
permModels  # Overall P-values for all K

# Save output to 
save(permModels, file="04.TCGA_Breast_Normal/Files/2016-08-04_Permutation_1000_age_BRCAnormal.RData")

# Separately by K
sapply(permModels, individual.p.values)

# Over just 1 to 5
sapply(permModels, summary, functional=make.minp.functional.subclasses(1:6))

# Over just 1 to 2
sapply(permModels, summary, functional=make.minp.functional.subclasses(1:2))

# Over just class 1 (should match above "separately by K" for K=1
sapply(permModels, summary, functional=make.minp.functional.subclasses(1))
sapply(permModels, summary, functional=make.minp.functional.subclasses(2))
sapply(permModels, summary, functional=make.minp.functional.subclasses(3))
sapply(permModels, summary, functional=make.minp.functional.subclasses(4))
sapply(permModels, summary, functional=make.minp.functional.subclasses(5))
sapply(permModels, summary, functional=make.minp.functional.subclasses(6))
sapply(permModels, summary, functional=make.minp.functional.subclasses(7))
sapply(permModels, summary, functional=make.minp.functional.subclasses(8))
sapply(permModels, summary, functional=make.minp.functional.subclasses(9))

# Plot for Permutation P-values
png("04.TCGA_Breast_Normal/Figures/TCGA_BRCA[n]_Permutation_Khat.png", height=7*300, width=7*300, res=300)
plot(x=c(2,3,4,5,6,7,8,9,10), y=-log10(sapply(permModels, individual.p.values)[ , 1]), pch=16, bty="n", col='red',
     ylab='-log10(P-value)', xlab='K (putative cell-types)', ylim=c(0,4.5))
lines(x=c(2,3,4,5,6,7,8,9,10), y=-log10(sapply(permModels, individual.p.values)[ , 1]), col='red', lwd = 2)
text(7, -log10(0.025), "Permutation P = 0.027", col='red', adj=c(0,-3), cex=1)
abline(h=-log10(0.05), col = "lightgray", lty = 2, lwd = 2)
text(7.25, -log10(0.04), "P < 0.05 threshold", col="lightgray", cex=1)
legend(2, 4, c("Age"), cex=1, lwd = 3, lty = 1, col=c('red'), bty='n')
dev.off()


