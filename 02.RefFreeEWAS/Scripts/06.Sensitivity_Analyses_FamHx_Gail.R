#####################################################
# Risk Factor-Related DNA methylation in Normal Breast Tissues 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 6: Sensitivity analyses using limma models testing association with Family Hx and Gail Score
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Packages used
library(data.table)
library(limma)
library(qvalue)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load Betas that have been processed by minfi pipeline
Komen_Betas <- data.frame(fread("01.Preprocessing/Data/Komen_BMIQ_Betas_17May2016.csv"), row.names=1)
Betas = Komen_Betas[ ,order(colnames(Komen_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("01.Preprocessing/Files/Komen_Manifest_20July2015.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]

# Inspect whether or not the samples are in the same order
all(names(Betas)==covars$Match_IDs)

# Keep all samples in analysis 
covars_Komen <- covars
rownames(covars_Komen) <- covars_Komen$Match_IDs
Betas_Komen <- data.matrix(Betas)
# Discount double check
all(names(Betas_Komen)==covars_Komen$Match_IDs)

# Cellular proportions estimated from RefFreeEWAS NMF. Optimal K (i.e., K-hat) = 6.
load("02.RefFreeEWAS/Files/Komen_Omega_BMIQ.RData")

# Same orientation for cell proportions and DNAme?
all(rownames(Komen_Omega)==names(Betas_Komen))

# Remove 1 cell-type to prevent multi-collinearity
Komen_om <- Komen_Omega[ , 1:5]
###########################################################
# Re-run CpG analysis with FamHx (remove 10 subjects with NAs)
###########################################################
# Sensitivity analysis for Family History and DNAmeth
samples_with_missing <- which(is.na(covars_Komen[ , 'BLDCAN']))
Komen_om_fam <- Komen_Omega[ -samples_with_missing, 1:5]

# Create a design matrix for covariates of interest
covars_Komen_fam <- covars_Komen[ -samples_with_missing, ]
XX_fam <- model.matrix(~Age+BMI+BLDCAN, data= covars_Komen_fam)
rownames(XX_fam) <- covars_Komen_fam$Match_IDs

# Convert beta-values to M-values for gaussian consideration
Betas_KomenM <- ifelse(Betas_Komen>=1, 1-1E-6, ifelse(Betas_Komen<=0, 1E-6, Betas_Komen))
Betas_KomenM <- log(Betas_KomenM)-log(1-Betas_KomenM)
colnames(Betas_KomenM) <- substring(colnames(Betas_KomenM), 2, nchar(colnames(Betas_KomenM))[1])

# Convert beta-values to M-values for gaussian consideration
Betas_KomenM_fam <- Betas_KomenM[ ,-samples_with_missing]
all(colnames(Betas_KomenM_fam)==rownames(XX_fam))

# Apply limma in unadjusted (Null) and cellular adjusted models
lf_Null_fam <-  eBayes(lmFit(Betas_KomenM_fam, XX_fam))
lf_Omega_fam <- eBayes(lmFit(Betas_KomenM_fam, cbind(XX_fam, Komen_om_fam)))

# Visualize the p-values
hist(lf_Omega_fam$p.value[ , 'BLDCANYes']) # uniform distribution of P-values

# Adjustment for multiple comparisons (note: use of q-value with threshold of 0.01)
q_BLDCANYes <- qvalue(lf_Omega_fam$p.value[ , 'BLDCANYes'], fdr.level = 0.01)
table(q_BLDCANYes$significant) # 0 Significant CpGs

# Calculate the difference in beta-coefficients between the unadj. and adj. model
Delta_Age <- lf_Null_fam$coef[,"BLDCANYes"] - lf_Omega_fam$coef[,"BLDCANYes"]

# Create data frames to feed to ggplot
RefFree_unadj_Data_Frame <- as.data.frame(cbind(lf_Null_fam$p.value, lf_Null_fam$coef, Delta_Age))
names(RefFree_unadj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "FamHx_Pval", "Inter2", "Age_coef", "BMI_coef", "FamHx_coef", "Delta_FamHx")
RefFree_adj_Data_Frame <- as.data.frame(cbind(lf_Omega_fam$p.value[, 1:4], lf_Omega_fam$coef[ ,1:4], Delta_Age))
names(RefFree_adj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "FamHx_Pval", "Inter2", "Age_coef", "BMI_coef", "FamHx_coef", "Delta_FamHx")

# ggplot for unadjusted model
unadj_limma = ggplot(RefFree_unadj_Data_Frame, aes(FamHx_coef, -log10(FamHx_Pval))) +
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_FamHx)) + 
  scale_color_gradient2(low = "grey", mid = "grey", high = "grey") + 
  labs(list(x = "FamHx coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Unadjusted"),
            color= "Delta\nlimma\nCoefficient")) + 
  xlim(-1.75, 1.75) + ylim(0, 20) + theme_bw()


# ggplot for adjusted model
adj_limma <- ggplot(RefFree_adj_Data_Frame, aes(FamHx_coef, -log10(FamHx_Pval))) + 
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_FamHx)) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  labs(list(x = "FamHx coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Adjusted"), 
            color = "Delta\nlimma\nCoefficient")) + 
  xlim(-1.75, 1.75) + ylim(0, 20) + theme_bw()


# Generate high-resolution heat maps as .pngs
library(gridExtra)
png("02.RefFreeEWAS/Figures/Figure1-Supplement-FamHX.png", height=7*300, width=7*300, res=300)
grid.arrange(unadj_limma, adj_limma, ncol = 2, nrow = 1)
dev.off()











###########################################################
# Re-run CpG analysis with FamHx (remove 55 subjects with NAs, less than 35 years old)
###########################################################
# Sensitivity analysis for Family History and DNAmeth
samples_with_missing_gail <- which(is.na(covars_Komen[ , 'Gail_Score']))
Komen_om_gail <- Komen_Omega[ -samples_with_missing_gail, 1:5]

# Create a design matrix for covariates of interest
covars_Komen_gail <- covars_Komen[ -samples_with_missing_gail, ]
XX_gail <- model.matrix(~Gail_Score, data = covars_Komen_gail)
rownames(XX_gail) <- covars_Komen_gail$Match_IDs

# Convert beta-values to M-values for gaussian consideration
Betas_KomenM <- ifelse(Betas_Komen>=1, 1-1E-6, ifelse(Betas_Komen<=0, 1E-6, Betas_Komen))
Betas_KomenM <- log(Betas_KomenM)-log(1-Betas_KomenM)
colnames(Betas_KomenM) <- substring(colnames(Betas_KomenM), 2, nchar(colnames(Betas_KomenM))[1])

# Convert beta-values to M-values for gaussian consideration
Betas_KomenM_gail <- Betas_KomenM[ ,-samples_with_missing_gail]
all(colnames(Betas_KomenM_gail)==rownames(XX_gail))

# Apply limma in unadjusted (Null) and cellular adjusted models
lf_Null_gail <-  eBayes(lmFit(Betas_KomenM_gail, XX_gail))
lf_Omega_gail <- eBayes(lmFit(Betas_KomenM_gail, cbind(XX_gail, Komen_om_gail)))

# Visualize the p-values
hist(lf_Omega_gail$p.value[ , 'Gail_Score']) # uniform distribution of P-values

# Adjustment for multiple comparisons (note: use of q-value with threshold of 0.01)
q_Gail <- qvalue(lf_Omega_gail$p.value[ , 'Gail_Score'], fdr.level = 0.01)
table(q_Gail$significant) # 0 Significant CpGs


