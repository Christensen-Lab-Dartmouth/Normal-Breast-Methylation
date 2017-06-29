###########################
#Genome-wide DNA methylation and its association with breast cancer risk factors
#Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate volcano plots for Limma unadjusted and adjusted models
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Packages used:
library(data.table)
library(limma)
library(qvalue)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load Betas that have been processed by minfi pipeline
Komen_Betas <- data.frame(fread("01.Preprocessing/Data/Komen_BMIQ_Betas_17May2016.csv"), row.names=1)
colnames(Komen_Betas) <- substring(colnames(Komen_Betas), 2, length(colnames(Komen_Betas)))
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
##############################################
# Limma approach for CpG-specific associations
##############################################
# Remove 1 cell-type to prevent multi-collinearity
Komen_om <- Komen_Omega[ , 1:5]

# Create a design matrix for covariates of interest
# covars_Komen$BMI_group <- ifelse(covars_Komen$BMI>=25.0, "overweight", "normal") 
XX <- model.matrix(~ Age + BMI + Pregnant, data= covars_Komen)
rownames(XX) <- covars_Komen$Match_IDs

# Convert beta-values to M-values for gaussian consideration
Betas_KomenM <- ifelse(Betas_Komen>=1, 1-1E-6, ifelse(Betas_Komen<=0, 1E-6, Betas_Komen))
Betas_KomenM <- log(Betas_KomenM)-log(1-Betas_KomenM)
all(colnames(Betas_KomenM)==rownames(XX))

# Apply limma in unadjusted (Null) and cellular adjusted models
lf_Null <-  eBayes(lmFit(Betas_KomenM, XX))
lf_Omega <- eBayes(lmFit(Betas_KomenM, cbind(XX, Komen_om)))

# Adjustment for multiple comparisons (note: use of q-value with threshold of 0.01)
q_age <- qvalue(lf_Omega$p.value[ , 'Age'], fdr.level = 0.01)
table(q_age$significant) # 787

# Save output results from EWAS with age
Age_DMPs <- names(which(q_age$qvalues<=0.01))
qThresh = max(q_age$pvalues[q_age$qvalues<=0.01])
qThresh  # 2.394087e-05 
multiple_correction <- -log10(qThresh)

# Create a data.frame for ggplot2
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")
annot$UCSC_RefGene_Name <- as.character(annot$UCSC_RefGene_Name)
annot$UCSC_RefGene_Group <- as.character(annot$UCSC_RefGene_Group)

# Extracted CpG and Gene region information
annot$Gene <- unlist(lapply(strsplit(annot$UCSC_RefGene_Name, ";"), '[', 1))
annot$Region <- unlist(lapply(strsplit(annot$UCSC_RefGene_Group, ";"), '[', 1))
annot$Gene_Region <- paste(annot$Gene, annot$Region, sep=":")
annot_trim <- annot[ , c('Name','Gene', 'Region', 'Gene_Region')]

# Generate Supplemental Table 2: annotation of 787 CpG sites
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
annot_supplemental = annot[Age_DMPs, c('CHR','MAPINFO','Gene', 'Region', 'Relation_to_UCSC_CpG_Island') ]
levels(annot_supplemental$Relation_to_UCSC_CpG_Island)[levels(annot_supplemental$Relation_to_UCSC_CpG_Island)==""] <- "OpenSea"
annot_supplemental$Relation_to_UCSC_CpG_Island <- gsub("^[N|S]_","", annot_supplemental$Relation_to_UCSC_CpG_Island) 
limma_celltype_adj = cbind(lf_Omega$p.value [, 2:4], lf_Omega$coefficients[ , 2:4])
annotated_results = merge(annot_supplemental, limma_celltype_adj, by="row.names")
write.csv(annotated_results, file="05.Genomic_Enrichment/Files/Age-Related-CpGs-Annotated.csv")

# Calculate the difference in beta-coefficients between the unadj. and adj. model
Delta_Age <- lf_Null$coef[,"Age"] - lf_Omega$coef[,"Age"]
RefFree_limma_adj_Komen <- as.data.frame(cbind(lf_Omega$p.value[ , 1:2], lf_Omega$coefficients[ ,1:2], Delta_Age))
names(RefFree_limma_adj_Komen) <- c("Inter1","Age_Pval","Inter2", "Age_coef", "DeltaAge")
Annotated_df <- merge(RefFree_limma_adj_Komen, annot_trim, by='row.names')
rownames(Annotated_df) <- Annotated_df$Row.names
Annotated_df$Row.names <- NULL

# Select the results that meet the threshold for significance
results = mutate(Annotated_df, sig=ifelse(Annotated_df$Age_Pval<2.394087e-05, "Q-value < 0.01", "Not Sig"))

# Figure 1. Genes with labeled regions
p = ggplot(results, aes(Age_coef, -log10(Age_Pval))) +
  geom_point(aes(col=sig)) + xlim(-.075, 0.075) + 
  xlab("Age coefficient - limma") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red"))
png("02.RefFreeEWAS/Figures/Figure1.png", height=7*300, width=7*300, pointsize = 12, res=300)
p + theme_bw() + geom_text_repel(data=filter(results, Age_Pval<1.3606e-17), aes(label=Gene_Region))
dev.off()

########################
# Delta Plot
########################
# Create data frames to feed to ggplot
RefFree_unadj_Data_Frame <- as.data.frame(cbind(lf_Null$p.value, lf_Null$coef, Delta_Age))
names(RefFree_unadj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "Parity_Pval", "Inter2", "Age_coef", "BMI_coef", "Parity_coef", "Delta_Age")
RefFree_adj_Data_Frame <- as.data.frame(cbind(lf_Omega$p.value[, 1:4], lf_Omega$coef[ ,1:4], Delta_Age))
names(RefFree_adj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "Parity_Pval", "Inter2", "Age_coef", "BMI_coef", "Parity_coef", "Delta_Age")

# ggplot for unadjusted model
unadj_limma = ggplot(RefFree_unadj_Data_Frame, aes(Age_coef, -log10(Age_Pval))) +
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_Age)) + 
  scale_color_gradient2(low = "grey", mid = "grey", high = "grey") + 
  labs(list(x = "Age coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Unadjusted"),
            color= "Delta\nlimma\nCoefficient")) + 
  xlim(-0.075, 0.075) + ylim(0, 20) + theme_bw()


# ggplot for adjusted model
adj_limma <- ggplot(RefFree_adj_Data_Frame, aes(Age_coef, -log10(Age_Pval))) + 
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_Age)) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  labs(list(x = "Age coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Adjusted"), 
            color = "Delta\nlimma\nCoefficient")) + 
  xlim(-0.075, 0.075) + ylim(0, 20) + theme_bw()

# Generate high-resolution heat maps as .pngs
library(gridExtra)
png("02.RefFreeEWAS/Figures/Figure1-Supplement-Age.png", height=7*300, width=7*300, res=300)
grid.arrange(unadj_limma, adj_limma, ncol = 2, nrow = 1)
dev.off()

############
#    BMI
############
Delta_BMI <- lf_Null$coef[,"BMI"] - lf_Omega$coef[,"BMI"]
# Create data frames to feed to ggplot
RefFree_unadj_Data_Frame <- as.data.frame(cbind(lf_Null$p.value, lf_Null$coef, Delta_BMI))
names(RefFree_unadj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "Parity_Pval", "Inter2", "Age_coef", "BMI_coef", "Parity_coef", "Delta_BMI")
RefFree_adj_Data_Frame <- as.data.frame(cbind(lf_Omega$p.value[, 1:4], lf_Omega$coef[ ,1:4], Delta_BMI))
names(RefFree_adj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "Parity_Pval", "Inter2", "Age_coef", "BMI_coef", "Parity_coef", "Delta_BMI")

# ggplot for unadjusted model
un_BMI <- ggplot(RefFree_unadj_Data_Frame, aes(BMI_coef, -log10(BMI_Pval))) + 
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_BMI)) + 
  scale_color_gradient2(low = "grey", mid = "grey", high = "grey") + 
  labs(list(x = "BMI coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Unadjusted"),
            color= "Delta\nlimma\nCoefficient")) + 
  xlim((-1 * 0.075), 0.075) + ylim(0, 20) + theme_bw()

# ggplot for adjusted model
ad_BMI <- ggplot(RefFree_adj_Data_Frame, aes(BMI_coef, -log10(BMI_Pval))) + 
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_BMI)) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  labs(list(x = "BMI coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Adjusted"), 
            color = "Delta\nlimma\nCoefficient")) + 
    xlim((-1 * 0.075), 0.075) + ylim(0, 20) + theme_bw()

# Generate high-resolution heat maps as .pngs
png("02.RefFreeEWAS/Figures/BMI_supplement.png", height=7*300, width=7*300, res=300)
grid.arrange(un_BMI, ad_BMI, ncol = 2, nrow = 1)
dev.off()

############
#    Parity
############
Delta_Parity <- lf_Null$coef[,"PregnantYes"] - lf_Omega$coef[,"PregnantYes"]
# Create data frames to feed to ggplot
RefFree_unadj_Data_Frame <- as.data.frame(cbind(lf_Null$p.value, lf_Null$coef, Delta_Parity))
names(RefFree_unadj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "Par_Pval", "Inter2", "Age_coef", "BMI_coef", "Par_coef", "Delta_Parity")
RefFree_adj_Data_Frame <- as.data.frame(cbind(lf_Omega$p.value[, 1:4], lf_Omega$coef[ ,1:4], Delta))
names(RefFree_adj_Data_Frame) <- c("Inter1","Age_Pval", "BMI_Pval", "Par_Pval", "Inter2", "Age_coef", "BMI_coef", "Par_coef", "Delta_Parity")

# ggplot for unadjusted model
un_Parity <- ggplot(RefFree_unadj_Data_Frame, aes(Par_coef, -log10(Par_Pval))) + 
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_Parity)) + 
  scale_color_gradient2(low = "grey", mid = "grey", high = "grey") + 
  labs(list(x = "Parity coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Unadjusted"),
            color= "Delta\nlimma\nCoefficient")) + 
    xlim(-1.75, 1.75) + ylim(0, 20) + theme_bw()

# ggplot for adjusted model
ad_Parity <- ggplot(RefFree_adj_Data_Frame, aes(Par_coef, -log10(Par_Pval))) + 
  geom_point(aes(colour = RefFree_adj_Data_Frame$Delta_Parity)) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  labs(list(x = "Parity coefficient - limma", 
            y = "-log10(P-value)", 
            title = paste("Cell Mixture\n", "Adjusted"),
            color= "Delta\nlimma\nCoefficient")) + 
  xlim(-1.75, 1.75) + ylim(0, 20) + theme_bw()

# Generate high-resolution heat maps as .pngs
png("02.RefFreeEWAS/Figures/Parity_supplement.png", height=7*300, width=7*300, res=300)
grid.arrange(un_Parity, ad_Parity, ncol = 2, nrow = 1)
dev.off()



