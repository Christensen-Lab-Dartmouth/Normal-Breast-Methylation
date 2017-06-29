###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear models with and and without adjustment for Omega in TCGA adjacent normal samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# NECESSARY libraries
library(RefFreeEWAS)
library(data.table)
library(limma)
library(ggplot2)
library(ggthemes)

# All age-related DMPs from Komen
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")

# Load Betas that have been processed by minfi
TCGA_Betas <- data.frame(fread("04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Betas_25July2016.csv"), row.names=1)
colnames(TCGA_Betas) <- substring(colnames(TCGA_Betas), 2, length(colnames(TCGA_Betas)))
Betas = TCGA_Betas[ ,order(colnames(TCGA_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("04.TCGA_Breast_Normal/Files/2016-08-04_TCGA_Normal_Breast_manifest.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]

# Discount double check
Betas_TCGA <- data.matrix(Betas[ , ])
all(names(Betas_TCGA) == covars$Match_IDs)
rownames(covars) <- covars$Match_IDs
XX <- model.matrix(~Age, data = covars)

# Convert Komen beta-values to M-values
Betas_TCGA_tmp <- ifelse(Betas_TCGA>=1,1-1E-6,ifelse(Betas_TCGA<=0,1E-6,Betas_TCGA))
Mvals_TCGA <- log(Betas_TCGA_tmp)-log(1-Betas_TCGA_tmp)

# Covariate files combined with Omega from RefFreeEWAS
TCGA_Omega <- read.csv("04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Cell_Proportions.csv", header=T)
# Renaming the putative cell types
colnames(TCGA_Omega)[11:17] <- c("Cell_type1", "Cell_type2","Cell_type3", "Cell_type4","Cell_type5", "Cell_type6",
                                 "Cell_type7")
# Fit limmas to DNA methylation data
all(colnames(Mvals_TCGA)==rownames(XX))
TCGA_lf_Null <-  eBayes(lmFit(Mvals_TCGA, XX))
TCGA_lf_Omega <- eBayes(lmFit(Mvals_TCGA, cbind(XX, (TCGA_Omega)[ , 11:16])))

# Not all CpGs will have passed QC in the TCGA data set (potentially probes were tossed out)
DMP_overlap <- Age_DMPs[which(Age_DMPs%in%rownames(TCGA_lf_Omega$p.value[ , ]))]

# Test the extent to which CpGs demonstrate an assoc with Age in TCGA
sum(TCGA_lf_Null$p.value[DMP_overlap, 'Age'] < 0.05) # 503 CpGs / 786 CpGs
sum(TCGA_lf_Omega$p.value[DMP_overlap , 'Age'] < 0.05) # 548 / 786 --- 70 % validated at nominal threshold!

# One of the hyper-methylated with age CpGs was missing in the TCGA data set
hyper_overlap <- hyper_Age_DMPs[which(hyper_Age_DMPs%in%rownames(TCGA_lf_Null$p.value[ , ]))]

# Determine whether hyper or hypo methylated CGs had greater reproducibility
sum(TCGA_lf_Omega$p.value[hyper_overlap, 'Age'] < 0.05 & TCGA_lf_Omega$coef[hyper_overlap,'Age'] > 0) # 410/544 CpGs. 75%
sum(TCGA_lf_Omega$p.value[hypo_Age_DMPs, 'Age'] < 0.05 & TCGA_lf_Omega$coef[hypo_Age_DMPs,'Age'] < 0) # 138/242 CpGs. 57%
sum(TCGA_lf_Omega$p.value[hyper_overlap, 'Age'] < (0.05/787) & TCGA_lf_Omega$coef[hyper_overlap,'Age'] > 0) # 313/544 CpGs. 57.5%
sum(TCGA_lf_Omega$p.value[hypo_Age_DMPs, 'Age'] < (0.05/787) & TCGA_lf_Omega$coef[hypo_Age_DMPs,'Age'] < 0) # 32/242 CpGs. 13.2%

######################################
#  Volcano Plots TCGA
######################################
RefFree_adj_TCGA_Data_Frame <- as.data.frame(cbind(TCGA_lf_Omega$p.value[DMP_overlap, 1:2], TCGA_lf_Omega$coefficients[DMP_overlap,1:2]))
names(RefFree_adj_TCGA_Data_Frame) <- c("Inter1","Age_Pval","Inter2", "Age_coef")

# Significance Threshold
# nom <- -log10(0.05)
# bonfer <- -log10(0.05/length(DMP_overlap))
# TCGA_limma = ggplot(RefFree_adj_TCGA_Data_Frame, aes(Age_coef, -log10(Age_Pval))) +
#  geom_point() 
# png("04.TCGA_Breast_Normal/Figures/2016-08-05_RefFree_TCGA_Age_volcano.png")
# TCGA_limma + theme_classic(base_size = 16) + 
# theme(axis.line.x = element_line(color="black"),
#      axis.line.y = element_line(color="black")) +
# geom_hline(yintercept = nom, color = "black", linetype = "dashed", size = 1.25) +
# geom_hline(yintercept = bonfer, color = "red", linetype = "dashed", size = 1.25) +
# xlim(-abs(max(RefFree_adj_TCGA_Data_Frame$Age_coef)), abs(max(RefFree_adj_TCGA_Data_Frame$Age_coef))) +
# labs(list(x = "M-value:subject age limma coefficient", 
#            y = "-log10(P-value)"))
# dev.off()

# Select the results that meet the threshold for significance
results = mutate(RefFree_adj_TCGA_Data_Frame, sig=ifelse(RefFree_adj_TCGA_Data_Frame$Age_Pval<0.05, "P-value < 0.05", "Not Sig"))

# Figure 1. Genes with labeled regions
p = ggplot(results, aes(Age_coef, -log10(Age_Pval))) +
  geom_point(aes(col=sig)) + xlim(-.075, 0.075) + 
  xlab("Age coefficient - limma") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red"))
png("04.TCGA_Breast_Normal/Figures/RefFree_TCGA_Age_volcano.png", height=7*300, width=7*300, pointsize = 12, res=300)
p + theme_bw() + geom_hline(yintercept=-log10(6.4E-05), lty="dashed")
# + geom_text_repel(data=filter(results, Age_Pval<1.25e-13), aes(label=Gene_Region))
dev.off()

######################################
#  Visualize top-hit for most highly significant in Komen population
######################################
TCGA_top_hit <- as.data.frame(cbind(Mvals_TCGA["cg07303143", ], covars$Age))
colnames(TCGA_top_hit) <- c("cg07303143", "Age")

# Plot the lowest P-value for associations with subject age
png("04.TCGA_Breast_Normal/Figures/top_CpG_Mvalue_TCGA.png", height=7*300, width=7*300, res=300)
p <- ggplot(TCGA_top_hit, aes(x = Age, y = cg07303143)) + geom_point() 
p + geom_smooth(method='lm', se=TRUE) + theme_classic(base_size = 16) +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) +
  labs(list(x = "Subject Age", 
            y = paste(colnames(TCGA_top_hit)[1], "(M-value)", sep=" ")))
dev.off()




