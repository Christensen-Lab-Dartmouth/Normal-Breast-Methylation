###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Identification of CpGs with differential methylation in TCGA tumors compared with adjacent-normal
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(limma)
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load TCGA Breast betas processed by minfi
load("08.Age_DMPs_Cancer/Files/TCGA_BRCA_Betas.RData")
tcga_betas_brca <- Betas
rm(Betas)
covariates <- read.table("08.Age_DMPs_Cancer/Files/BRCAtarget_covariates.csv", header = T, sep = ",", stringsAsFactors = F)
# Subset covariate file
rownames(covariates) <- paste(covariates$Slide, covariates$Array, sep='_')
covariates <- covariates[!is.na(covariates$age.Dx), ]
Mets_indices <- which(covariates$sample.type=="Metastatic")
covariates <- covariates[-Mets_indices, ]
covariates$tissue_status <- rep(0, length(covariates$sample.type))
covariates$tissue_status[covariates$sample.type=="Primary Tumor"] <- 1

# Compare only TCGA adjacent normal and primary tumors (n=97)
covariates_matched <- covariates[intersect(rownames(covariates), colnames(tcga_betas_brca)),]
tcga_breast_matched <- tcga_betas_brca[ ,intersect(rownames(covariates), colnames(tcga_betas_brca))]
rm(tcga_betas_brca)
# Are the samples the same in both the TCGA betas and the covariates
all(colnames(tcga_breast_matched)==rownames(covariates_matched))

# Betas_TCGA <- Betas_TCGA[Age_DMPs[Age_DMPs%in%rownames(Betas_TCGA)], ]
# Create a design matrix for covariates of interest
table(covariates_matched$sample.type)
covariates_matched$sample.type <- factor(covariates_matched$sample.type, levels=c("Solid Tissue Normal", "Primary Tumor"))
covariates_matched$age.Dx <- as.numeric(covariates_matched$age.Dx)
XX <- model.matrix(~sample.type+age.Dx, data = covariates_matched)

# par(bty='l')
# boxplot(Betas_TCGA["cg18451114", XX[ , 'sample.typePrimary Tumor']==0], Betas_TCGA["cg18451114", XX[ , 'sample.typePrimary Tumor']==1], ylim=c(0,1), cex.axis = 1.5)
# points(x= jitter(rep(1,97), factor=3), Betas_TCGA["cg18451114", XX[ , 'sample.typePrimary Tumor']==0], col="red", pch=16)
# points(x= jitter(rep(2, 749), factor=1.5), Betas_TCGA["cg18451114", XX[ , 'sample.typePrimary Tumor']==1], col="red", pch=16)
# cg18451114 : 6.337548e-37

# Convert beta-values to M-values for gaussian consideration
tcga_breast_matched <- ifelse(tcga_breast_matched>=1,1-1E-6,ifelse(tcga_breast_matched<=0,1E-6,tcga_breast_matched))
Betas_TCGAM <- log(tcga_breast_matched)-log(1-tcga_breast_matched)
all(colnames(Betas_TCGAM)==rownames(XX))

# Apply limma in unadjusted (Null) and cellular adjusted models
lf_Null <-  eBayes(lmFit(Betas_TCGAM, XX))
save(lf_Null, file="/Users/kevinjohnson/Komen_Normal_5mC/08.Age_DMPs_Cancer/Files/TCGA_BRCA_limma_models.RData")

hist(lf_Null$p.value[ , 'sample.typePrimary Tumor'])
sum(lf_Null$p.value[ , 'sample.typePrimary Tumor'] < (0.05/dim(lf_Null$p.value)[1])) # 170,998 CpGs

# We are primarily interested in those age-related CpGs
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
overlap_CpGs <- Age_DMPs[which(Age_DMPs%in%rownames(lf_Null$p.value))]

sum(lf_Null$p.value[overlap_CpGs, 'sample.typePrimary Tumor'] < (0.05)) # 642/787

# Restrict the annotation file to only those CpGs used in the Komen analysis
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")
annot_overlap <- annot[overlap_CpGs, ] 
annot_results <- annot[rownames(lf_Null), ] 
levels(annot_overlap$Relation_to_UCSC_CpG_Island)[levels(annot_overlap$Relation_to_UCSC_CpG_Island)==""] <- "OpenSea"
levels(annot_results$Relation_to_UCSC_CpG_Island)[levels(annot_results$Relation_to_UCSC_CpG_Island)==""] <- "OpenSea"
# Collapse "North" and "South" nomenclature for CpG islands
annot_overlap$Relation_to_UCSC_CpG_Island <- gsub("^[N|S]_","", annot_overlap$Relation_to_UCSC_CpG_Island) 
annot_results$Relation_to_UCSC_CpG_Island <- gsub("^[N|S]_","", annot_results$Relation_to_UCSC_CpG_Island) 

# Randomly select CpGs with the same distribution of CpG island for ks.test
table(annot_overlap$Relation_to_UCSC_CpG_Island) # Island = 290, OpenSea = 224, Shelf = 57, Shore = 216
table(annot_results$Relation_to_UCSC_CpG_Island)
set.seed(451)
Island = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="Island", ]), 290)
set.seed(451)
OpenSea = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="OpenSea", ]), 224)
set.seed(451)
Shelf = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="Shelf", ]), 57)
set.seed(451)
Shore = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="Shore", ]), 216)
TCGA_random_CpGs = c(Island, OpenSea, Shelf, Shore)
# Prohibit randomly selected CpGs from being age-related
sum(TCGA_random_CpGs%in%overlap_CpGs) # none of the random CpGs are age-related


# Select the results that meet the threshold for significance
Limma_Results <- as.data.frame(cbind((lf_Null$coef[overlap_CpGs, 'sample.typePrimary Tumor']), lf_Null$p.value[overlap_CpGs, 'sample.typePrimary Tumor']))
results = mutate(Limma_Results, sig=ifelse(Limma_Results$V2<0.05, "P-value < 0.05", "Not Sig"))
colnames(results) = c("Coefficient", "P_value", "limma_model") 

# Figure for locus-specific TCGA vs. normal
p = ggplot(results, aes(Coefficient, -log10(P_value))) +
  geom_point(aes(col=limma_model))  + xlim(-2.5, 2.5) +
  xlab("IDC coefficient - limma") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red")) + 
  geom_hline(yintercept = -log10(.05), color = "red", linetype = "dashed", size = 1.3)
png("08.Age_DMPs_Cancer/Figures/TCGA_Normal_Volcano.png", height=7*250, width=7*300, pointsize = 16, res=300)
p + theme_bw() + labs(color="limma model \n significance") + theme(legend.key = element_blank()) + ggtitle("Tumor(n=749) -Normal(n=97) Differential \n DNA methylation (TCGA)")
dev.off()

# Sample to examine differences in coefficients for boxplot
random_sites <- lf_Null$coef[TCGA_random_CpGs, 'sample.typePrimary Tumor']
age_related_sites <- Limma_Results$V1
coef_compare <- cbind(age_related_sites, random_sites)

# kolmogorov-smirnov test between age-related and randomly selected CpGs
ks.test(age_related_sites, random_sites) # < 2.2E-16

# Provide supplemental plot for differences in coefficients
library(reshape)
coef_compare_long = melt(coef_compare, id=c("age_related_sites","random_sites"))
colnames(coef_compare_long) = c("CpG", "CpG_set", "limma_coefficient")
png("08.Age_DMPs_Cancer/Figures/TCGA_Normal_Coef_Density.png", height=7*300, width=7*300, res=300)
ggplot(coef_compare_long, aes(x = limma_coefficient, fill = CpG_set)) + geom_density(alpha = 0.5) + theme_bw() +
  xlab("limma coefficient") + ylab("Density")
dev.off()

# Apply kolmogorov-smirnov test to distribution of P-values
random_sites_pval <- -log10(lf_Null$p.value[TCGA_random_CpGs, 'sample.typePrimary Tumor'])
age_related_sites_pval <- -log10(Limma_Results$V2)

# kolmogorov-smirnov test to -log10()
ks.test(age_related_sites_pval, random_sites_pval) # P-value = 1.105e-13

# Plot distributions for supplemental figures
pval_compare <- cbind(age_related_sites_pval, random_sites_pval)
pval_compare_long = melt(pval_compare, id=c("age_related_sites","random_sites"))
colnames(pval_compare_long) = c("CpG", "CpG_set", "limma_pvalue")
png("08.Age_DMPs_Cancer/Figures/TCGA_Normal_Pval_Density.png", height=7*275, width=7*275, res=300)
ggplot(pval_compare_long, aes(x = limma_pvalue, fill = CpG_set)) + geom_density(alpha = 0.5) + theme_bw() +
  xlab("-log10(P-value)") + ylab("Density")
dev.off()

# Plot the cdf of the age-related and random P-values
age_pval_cdf = ecdf(age_related_sites_pval)
Random_pval_cdf = ecdf(random_sites_pval)
png("08.Age_DMPs_Cancer/Figures/TCGA_Normal_Pvalue_ecdf.png", height=7*300, width=7*275, pointsize = 12, res=300)
plot(age_pval_cdf, main="Tumor-Normal Differential \n DNA methylation (TCGA)", 
     xlab="-log10(P-value)", ylab="Cumulative proportion - P-value", col="black")
plot(Random_pval_cdf, add=TRUE, col="gray")
legend(40, 0.5, c("Komen age-related CpGs", "Random CpGs"), lwd=3, bty="n", col=c("black","gray"))
text(50, 0.3, "Kolmogorov-Smirnov test P = 1.1E-13")     
dev.off()


random_sites_coef <- lf_Null$coef[TCGA_random_CpGs, 'sample.typePrimary Tumor']
age_related_sites_coef <- Limma_Results$V1
ks.test(age_related_sites_coef, random_sites_coef) # < 2.2E-16

age_coef_cdf = ecdf(age_related_sites_coef)
random_coef_cdf = ecdf(random_sites_coef)
png("08.Age_DMPs_Cancer/Figures/TCGA_Normal_Coef_ecdf.png", height=7*300, width=7*300, pointsize = 12, res=300)
plot(age_coef_cdf, main="Tumor-Normal Differential \n DNA methylation (TCGA)", 
     xlab="limma coefficient", ylab="Cumulative Proportion - Coefficient", col="black")
plot(random_coef_cdf, add=TRUE, col="gray")
legend(.25, 0.5, c("Komen age-related CpGs", "Random CpGs"), lwd=3, bty="n", col=c("black","gray"))
text(1.5, 0.3, "Kolmogorov-Smirnov \n test P = <2.2E-16")     
dev.off()

