###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Volcano plot for adjusted and unadjusted RefFreeEWAS models in normal breast tissue (Komen)
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(limma)
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load Betas that have been processed by minfi
load("08.Age_DMPs_Cancer/Files/NHMN_DCIS_Betas_BMIQ.RData")
Betas <- Betas_BMIQ[ ,order(colnames(Betas_BMIQ), decreasing=T)]

NHMN_Betas <- data.frame(fread("08.Age_DMPs_Cancer/Files/DCIS_BMIQ_Betas.csv"), row.names=1)
colnames(NHMN_Betas) <- substring(colnames(NHMN_Betas), 2, length(colnames(NHMN_Betas)))
Betas = NHMN_Betas[ ,order(colnames(NHMN_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("08.Age_DMPs_Cancer/Files/DCIS_covariates.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T), ]

# Inspect whether or not the samples are in the same order
all(colnames(Betas)==covars$Match_IDs)

# Keep all samples in analysis 
covars_DCIS <- covars
Betas_DCIS <- data.matrix(Betas[, ])

# Discount double check
all(names(Betas_DCIS)==covars_DCIS$Match_IDs)

# Create a design matrix for covariates of interest
XX <- model.matrix(~Sample_Group + agepr, data= covars_DCIS)
rownames(XX) <- covars_DCIS$Match_IDs

# Convert beta-values to M-values for gaussian consideration
Betas_DCISM <- ifelse(Betas_DCIS>=1,1-1E-6,ifelse(Betas_DCIS<=0,1E-6,Betas_DCIS))
Betas_DCISM <- log(Betas_DCISM)-log(1-Betas_DCISM)
all(colnames(Betas_DCISM)==rownames(XX))

# Apply limma in unadjusted (Null) and cellular adjusted models
lf_Null <-  eBayes(lmFit(Betas_DCISM, XX))

hist(lf_Null$p.value[ , 'Sample_GroupDCIS'])
# Bonferroni significance threshold
sum(lf_Null$p.value[ , 'Sample_GroupDCIS'] < (0.05/dim(lf_Null$p.value)[1])) # 334 CpGs

# Overlap between significant CpGs and CpGs passing QC measures in the DCIS study
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
overlap_CpGs <- Age_DMPs[Age_DMPs%in%rownames(lf_Null$p.value)]
# At the nominal level
sum(lf_Null$p.value[overlap_CpGs, 'Sample_GroupDCIS'] < 0.05) # 268/775 CpGs measured

# Sample absolute value of coefficients between randomly sampled CpGs. Are age-related CpGs larger?

#wilcotestgenerator <- function(n) {
  # note that here we have a "age-related" group and a 
  # randomly selected group. The age-related remains the same (n=787)
#  age_related_sites <- abs(lf_Null$coefficients[overlap_CpGs , 'Sample_GroupDCIS'])
#  random_sites <- abs(sample(lf_Null$coefficients[ , 'Sample_GroupDCIS'], n))
#  wilco_stat <- wilcox.test(age_related_sites, random_sites) 
#  return(wilco_stat$p.value)
# }
# set.seed(123)
# wilco_tests_DCIS <- replicate(1000, wilcotestgenerator(775))
# summary(wilco_tests_DCIS) # median 6.14E-08


# Restrict the annotation file to only those CpGs used in the Komen analysis
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")
annot_overlap <- annot[overlap_CpGs, ] 
annot_results <- annot[rownames(Betas_DCISM), ] 
levels(annot_overlap$Relation_to_UCSC_CpG_Island)[levels(annot_overlap$Relation_to_UCSC_CpG_Island)==""] <- "OpenSea"
levels(annot_results$Relation_to_UCSC_CpG_Island)[levels(annot_results$Relation_to_UCSC_CpG_Island)==""] <- "OpenSea"
# Collapse "North" and "South" nomenclature for CpG islands
annot_overlap$Relation_to_UCSC_CpG_Island <- gsub("^[N|S]_","", annot_overlap$Relation_to_UCSC_CpG_Island) 
annot_results$Relation_to_UCSC_CpG_Island <- gsub("^[N|S]_","", annot_results$Relation_to_UCSC_CpG_Island) 

table(annot_overlap$Relation_to_UCSC_CpG_Island) # Island = 289, OpenSea = 218, Shelf = 54, Shore = 214
table(annot_results$Relation_to_UCSC_CpG_Island)
set.seed(451)
Island = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="Island", ]), 289)
set.seed(451)
OpenSea = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="OpenSea", ]), 218)
set.seed(451)
Shelf = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="Shelf", ]), 54)
set.seed(451)
Shore = sample(rownames(annot_results[annot_results$Relation_to_UCSC_CpG_Island=="Shore", ]), 214)
random_CpGs = c(Island, OpenSea, Shelf, Shore)
sum(random_CpGs%in%overlap_CpGs) # none of the random CpGs are age-related

# Provide visualizations for 
# Select the results that meet the threshold for significance
Limma_Results <- as.data.frame(cbind((lf_Null$coef[overlap_CpGs, 'Sample_GroupDCIS']), lf_Null$p.value[overlap_CpGs, 'Sample_GroupDCIS']))
results = mutate(Limma_Results, sig=ifelse(Limma_Results$V2<0.05, "P-value < 0.05", "Not Sig"))
colnames(results) = c("Coefficient", "P_value", "limma_model") 
# Figure for locus-specific DCIS vs. normal
p = ggplot(results, aes(Coefficient, -log10(P_value))) +
  geom_point(aes(col=limma_model))  + xlim(-2.5, 2.5) +
  xlab("DCIS coefficient - limma") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red")) + 
  geom_hline(yintercept = -log10(.05), color = "red", linetype = "dashed", size = 1.3)
png("08.Age_DMPs_Cancer/Figures/DCIS_Normal_Volcano.png", height=7*250, width=7*300, pointsize = 12, res=300)
p + theme_bw() + labs(color="limma model \n significance") + theme(legend.key = element_blank()) + ggtitle("Tumor(n=40) -Normal(n=15) Differential \n DNA methylation (DCIS)")
dev.off()

# Compare the P-values from the randomly selected set of loci (above) to the age-related CpGs
random_sites_pval <- -log10(lf_Null$p.value[random_CpGs, 'Sample_GroupDCIS'])
age_related_sites_pval <- -log10(Limma_Results$V2)
ks.test(age_related_sites_pval, random_sites_pval)
pval_compare <- cbind(age_related_sites_pval, random_sites_pval)

# Plot the cdf of the age-related and random P-values
age_pval_cdf = ecdf(age_related_sites_pval)
Random_pval_cdf = ecdf(random_sites_pval)
png("08.Age_DMPs_Cancer/Figures/DCIS_Normal_Pvalue_ecdf.png", height=7*300, width=7*275, pointsize = 12, res=300)
plot(age_pval_cdf, main="Tumor-Normal Differential \n DNA methylation (DCIS)", 
     xlab="-log10(P-value)", ylab="Cumulative proportion - P-value", col="black")
plot(Random_pval_cdf, add=TRUE, col="gray")
legend(4, 0.5, c("Komen age-related CpGs", "Random CpGs"), lwd=3, bty="n", col=c("black","gray"))
text(5, 0.3, "Kolmogorov-Smirnov test P = 3.0E-03")     
dev.off()

# Compare the limma coefficients from the set of randomly selected CpGs and age-related (Supplementary Figure)
random_sites_coef <- lf_Null$coef[random_CpGs, 'Sample_GroupDCIS']
age_related_sites_coef <- Limma_Results$V1
ks.test(age_related_sites_coef, random_sites_coef)

age_coef_cdf = ecdf(age_related_sites_coef)
random_coef_cdf = ecdf(random_sites_coef)
png("08.Age_DMPs_Cancer/Figures/DCIS_Normal_Coef_ecdf.png", height=7*300, width=7*300, pointsize = 12, res=300)
plot(age_coef_cdf, main="Tumor-Normal Differential \n DNA methylation (DCIS)", 
     xlab="limma coefficient", ylab="Cumulative Proportion - Coefficient", col="black")
plot(random_coef_cdf, add=TRUE, col="gray")
legend(.5, 0.5, c("Komen age-related CpGs", "Random CpGs"), lwd=3, bty="n", col=c("black","gray"))
text(1.5, 0.3, "Kolmogorov-Smirnov \n test P = 8.9E-15")     
dev.off()


# Plot distribution of P-values
pval_compare_long = melt(pval_compare, id=c("age_related_sites","random_sites"))
colnames(pval_compare_long) = c("CpG", "CpG_set", "limma_pvalue")
png("08.Age_DMPs_Cancer/Figures/DCIS_Normal_Pval_Density.png", height=7*300, width=7*300, res=300)
ggplot(pval_compare_long, aes(x = limma_pvalue, fill = CpG_set)) + geom_density(alpha = 0.5) + theme_bw() +
  xlab("-log10(P-value)") + ylab("Density")
dev.off()

# Plot the distribution of coefficients from age-related and the randomly selected CpGs
random_sites <- lf_Null$coef[random_CpGs, 'Sample_GroupDCIS']
age_related_sites <- Limma_Results$V1
coef_compare <- cbind(age_related_sites, random_sites)
ks.test(age_related_sites, random_sites)

library(reshape)
coef_compare_long = melt(coef_compare, id=c("age_related_sites","random_sites"))
colnames(coef_compare_long) = c("CpG", "CpG_set", "limma_coefficient")
png("08.Age_DMPs_Cancer/Figures/DCIS_Normal_Coef_Density.png", height=7*300, width=7*300, res=300)
ggplot(coef_compare_long, aes(x = limma_coefficient, fill = CpG_set)) + geom_density(alpha = 0.5) + theme_bw() +
  xlab("limma coefficient") + ylab("Density")
dev.off()






