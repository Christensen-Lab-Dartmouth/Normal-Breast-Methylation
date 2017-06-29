###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear models with and and without adjustment for Omega in NDRI samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(limma)
library(data.table)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC")

# All age-related DMPs from Komen
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")

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

# Examine all CpGs
NDRI_Betas <- as.matrix(Betas[, ]) 

# Discount double check
all(names(NDRI_Betas)== covars$Match_IDs)
rownames(covars) <- covars$Match_IDs
XX <- model.matrix(~Sample_Age + Sample_BMI, data= covars)

# Convert Komen beta-values to M-values
Betas_NDRIM <- ifelse(NDRI_Betas>=1,1-1E-6,ifelse(NDRI_Betas<=0,1E-6,NDRI_Betas))
Betas_NDRIM <- log(Betas_NDRIM)-log(1-Betas_NDRIM)

# Covariate files combined with Omega from RefFreeEWAS
NDRI_Omega <- read.csv("03.NDRI_Validation/Files/NDRI_BMIQ_Cell_Proportions.csv", header=T)
colnames(NDRI_Omega)[10:11] <- c("Cell_type1", "Cell_type2")

all(colnames(Betas_NDRIM)==rownames(XX))
lf_Null_NDRI <-  eBayes(lmFit(Betas_NDRIM, XX))
lf_Omega_NDRI <- eBayes(lmFit(Betas_NDRIM, cbind(XX, NDRI_Omega[,10])))

# Bonferonni threshold for cell-type adj and unadj models
sum(lf_Null_NDRI$p.value[ , 'Sample_Age'] < (0.05/dim(lf_Null_NDRI$p.value)[1])) # 75 
sum(lf_Omega_NDRI$p.value[ , 'Sample_Age'] < (0.05/787)) # 578 


# Overlap of significance at nominal-threshold...81.5% overlap with directionally consistent relationships
sum(lf_Omega_NDRI$p.value[hyper_Age_DMPs, 'Sample_Age'] < 0.05 & lf_Omega_NDRI$coef[hyper_Age_DMPs,'Sample_Age'] > 0) # 372 / 545 CpGs
sum(lf_Omega_NDRI$p.value[hypo_Age_DMPs, 'Sample_Age'] < 0.05 & lf_Omega_NDRI$coef[hypo_Age_DMPs,'Sample_Age'] < 0) # 17 / 242 CpGs
sum(lf_Omega_NDRI$p.value[hyper_Age_DMPs, 'Sample_Age'] < (0.05/787) & lf_Omega_NDRI$coef[hyper_Age_DMPs,'Sample_Age'] > 0) # 372 / 545 CpGs
sum(lf_Omega_NDRI$p.value[hypo_Age_DMPs, 'Sample_Age'] < (0.05/787) & lf_Omega_NDRI$coef[hypo_Age_DMPs,'Sample_Age'] < 0) # 17 / 242 CpGs

######################################
#  Volcano Plots NDRI
######################################
RefFree_adj_NDRI_Data_Frame <- as.data.frame(cbind(lf_Omega_NDRI$p.value[Age_DMPs, 1:2], lf_Omega_NDRI$coef[Age_DMPs, 1:2]))
names(RefFree_adj_NDRI_Data_Frame) <- c("Inter1","Age_Pval","Inter2", "Age_coef")

# Significance Threshold
# nom <- -log10(0.05)
# bonfer <- -log10(0.05/length(Age_DMPs))
# NDRI_limma = ggplot(RefFree_adj_NDRI_Data_Frame, aes(Age_coef, -log10(Age_Pval))) +
#  geom_point() 
# png("03.NDRI_Validation/Figures/2016-08-05_RefFree_NDRI_Age_volcano.png", height=7*300, width=7*300, res=300)
# NDRI_limma + theme_classic(base_size = 16) + 
#  theme(axis.line.x = element_line(color="black"),
#        axis.line.y = element_line(color="black")) +
#  geom_hline(yintercept = nom, color = "black", linetype = "dashed", size = 1.25) +
#  geom_hline(yintercept = bonfer, color = "red", linetype = "dashed", size = 1.25) +
#  xlim(-abs(max(RefFree_adj_NDRI_Data_Frame$Age_coef)), abs(max(RefFree_adj_NDRI_Data_Frame$Age_coef))) +
#  labs(list(x = "M-value:subject age limma coefficient", 
#            y = "-log10(P-value)"))
# dev.off()


# Select the results that meet the threshold for significance
results = mutate(RefFree_adj_NDRI_Data_Frame, sig=ifelse(RefFree_adj_NDRI_Data_Frame$Age_Pval<0.05, "P-value < 0.05", "Not Sig"))

# Figure 1. Genes with labeled regions
p = ggplot(results, aes(Age_coef, -log10(Age_Pval))) +
  geom_point(aes(col=sig)) + xlim(-.075, 0.075) + 
  xlab("Age coefficient - limma") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red"))
png("03.NDRI_Validation/Figures/RefFree_NDRI_Age_volcano.png", height=7*300, width=7*300, pointsize = 12, res=300)
p + theme_bw() + geom_hline(yintercept=-log10(6.4E-05), lty="dashed")
# + geom_text_repel(data=filter(results, Age_Pval<1.25e-13), aes(label=Gene_Region))
dev.off()


######################################
#  Visualize top-hit for most highly significant in Komen population
######################################
NDRI_top_hit <- as.data.frame(cbind(Betas_NDRIM["cg07303143", ], covars$Sample_Age))
colnames(NDRI_top_hit) <- c("cg07303143", "Age")

# Plot the lowest P-value for associations with subject age
png("03.NDRI_Validation/Figures/top_CpG_Mvalue_NDRI.png", height=7*300, width=7*300, pointsize = 12, res=300)
p <- ggplot(NDRI_top_hit, aes(x = Age, y = cg07303143)) + geom_point() 
p + geom_smooth(method='lm', se=TRUE) + theme_classic(base_size = 16) +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) +
  labs(list(x = "Subject Age", 
            y = paste(colnames(NDRI_top_hit)[1], "(M-value)", sep=" ")))
dev.off()


