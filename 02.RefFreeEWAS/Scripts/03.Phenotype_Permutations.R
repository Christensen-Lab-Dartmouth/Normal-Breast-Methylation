###########################
# Breast cancer risk factors are associated with DNA methylation in non-diseased breast tissue independent of cell type
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Permutation test for determining associations with metadata across classes
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Packages
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
Betas_Komen <- data.matrix(Betas)
# Samples are ordered the same in both matrices
all(names(Betas_Komen)== covars_Komen$Match_IDs)

# Build Design Matrix
DESIGNX <- model.matrix(~Age + BMI + Pregnant, data= covars_Komen)

# Load results from RefFree
load("02.RefFreeEWAS/Files/RefFree-2-18May2016.RData")
# Keq_BMIQ <- RefFreeCellMix(Betas_Komen, mu0=RefFreeCellMixInitialize(Betas_Komen, K=6))
# Komen_Omega <- Keq_BMIQ$Omega
# save(Komen_Omega, file="02.RefFreeEWAS/Files/Komen_Omega_BMIQ.RData")
load("02.RefFreeEWAS/Files/Komen_Omega_BMIQ.RData")

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
### EXAMPLES from Andy:
# BREAST TUMORS (#1) (27K)
# XBLOCK <- list(Tumor=2:4, Age=5, Size=6)
#
# BREAST TUMORS (#2) (27K)
# XBLOCK <- list(Tumor=2:6, Age=7)

# Normal Breast (450K)
XBLOCK <- list(Age=2, BMI=3, Pregnant=4)
XX <- model.matrix(~Age + BMI + Pregnant, data= covars_Komen)
Komen_om_List <- sapply(Komen_RefFree_Array, "[[", 2) 

# Inclusion of parity (ever pregnant? YES or NO)
permModels <- lapply(XBLOCK, function(u){
  cat(colnames(XX)[u],"\n")
  omegArrayPermutation(Komen_om_List[1:9], XX, R=1000, covariates=u)
})
permModels  # Overall P-values for all K

# save(permModels, file="02.RefFreeEWAS/Files/Permutation_1000_BMIQ_age_BMI_parity.RData")

# Separately by K
sapply(permModels, individual.p.values)

# Over just 1 to 5
sapply(permModels, summary, functional= make.minp.functional.subclasses(1:6))

# Over just 1 to 2
sapply(permModels, summary, functional= make.minp.functional.subclasses(1:2))

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

# Plot for Permutation P-values (Quick hack. Values are specific to my data set.)
ppi <- 300
png("02.RefFreeEWAS/Figures/Phenotype_Assoc_Permutation.png", height=7*ppi, width=7*ppi, res=ppi)
plot(x=c(2,3,4,5,6,7,8,9,10), y=-log10(c(0.203, 0.010, 0.004, 0.0001, 0.001, 0.0001, 0.001, 0.007, 0.001)), pch=16, bty="n", col='red',
     ylab='-log10(P-value)', xlab='K (putative cell-types)', ylim=c(0,4.5))
lines(x=c(2,3,4,5,6,7,8,9,10), y=-log10(c(0.203, 0.010, 0.004, 0.0001, 0.001, 0.0001, 0.001, 0.007, 0.001)), col='red')
points(x=c(2,3,4,5,6,7,8,9,10), y=-log10(sapply(permModels, individual.p.values)[ , 2]), pch=16, col='blue')
lines(x=c(2,3,4,5,6,7,8,9,10), y=-log10(sapply(permModels, individual.p.values)[ , 2]), col='blue')
points(x=c(2,3,4,5,6,7,8,9,10), y=-log10(sapply(permModels, individual.p.values)[ , 3]), pch=16, col='purple')
lines(x=c(2,3,4,5,6,7,8,9,10), y=-log10(sapply(permModels, individual.p.values)[ , 3]), col='purple')
text(5, -log10(0.002), "Perm P-val = 0.002", col='red', adj=c(0,-3), cex=0.8)
text(5, -log10(0.081), "Perm P-val = 0.08", col='blue', adj=c(0,-3.5), cex=0.8)
text(5, -log10(0.31), "Perm P-val = 0.31", col='purple', adj=c(0,-2), cex=0.8)
abline(h=-log10(0.05), col = "lightgray", lty=2)
text(9.5, -log10(0.04), "P<0.05 threshold", col="lightgray", cex=0.75)
legend(2, 4, c("Age","BMI", "Parity"), cex=1, lwd=2, lty = 1, col=c('red', 'blue', 'purple'), bty='n')
dev.off()

