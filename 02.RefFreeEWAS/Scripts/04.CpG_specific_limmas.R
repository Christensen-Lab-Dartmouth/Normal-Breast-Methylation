#####################################################
# Breast cancer risk factors are associated with DNA methylation in non-diseased breast tissue independent of cell type
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4: Limmas for determining CpG-specific associations for models unadjusted and adjusted for cell-type
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Packages used:
library(data.table)
library(limma)
library(qvalue)

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
library(corrgram)
corrgram(Komen_om, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="Cellular Proportions (K hat)")

# DIAGNOSTICS: Assess whether any two cell-types are highly correlated with one another. 
cor.test(Komen_om[ , '1'], Komen_om[ ,'4']) # PCC = -0.58 , P-val = 1.4E-10
cor.test(Komen_om[ , '4'], Komen_om[ ,'5']) # PCC = 0.31 , P-val = 0.0015 

ppi <- 300
png("02.RefFreeEWAS/Figures/cell_proportions_limma_model_boxplot.png", height=7*ppi, width=7*ppi, res=ppi)
par(bty='l')
boxplot(Komen_om, pch= 19, col='pink', xlab="Putative Cell-Types", ylab="Relative Proportions", main="Komen Normal Breast Tissues \n n=100")
dev.off()

# Is age correlated with BMI in this data set?
plot(covars_Komen$Age, covars_Komen$BMI, pch=16)
cor.test(covars_Komen$Age, covars_Komen$BMI) # PCC = 0.14, P-val = 0.16

# Create a design matrix for covariates of interest
# It is possible to treatment BMI as a binary variable with a threshold of 25, which is classified as ovrweight
# covars_Komen$BMI_group <- ifelse(covars_Komen$BMI>=25.0, "overweight", "normal") 
XX <- model.matrix(~ Age + BMI + Pregnant, data= covars_Komen)
rownames(XX) <- covars_Komen$Match_IDs

# Convert beta-values to M-values for gaussian considerations
Betas_KomenM <- ifelse(Betas_Komen>=1, 1-1E-6, ifelse(Betas_Komen<=0, 1E-6, Betas_Komen))
Betas_KomenM <- log(Betas_KomenM)-log(1-Betas_KomenM)
all(colnames(Betas_KomenM)==rownames(XX))

# Apply limma in unadjusted (Null) and cellular adjusted models
lf_Null <-  eBayes(lmFit(Betas_KomenM, XX))
lf_Omega <- eBayes(lmFit(Betas_KomenM, cbind(XX, Komen_om))) # adjusted-model

#####################################################
# Enumerate CpGs that exceed statistical significance
#####################################################
# Glance at top hit
which.min(lf_Omega$p.value[ , 'Age']) # cg07303143

Komen_top_hit <- as.data.frame(cbind(Betas_KomenM["cg06458239", ], covars_Komen$Age))
colnames(Komen_top_hit) <- c("cg07303143", "Age")

# Plot the lowest P-value for associations with subject age
library(ggplot2)
png("02.RefFreeEWAS/Figures/top_CpG_Mvalue_Komen.png",height=7*300, width=7*300, res=300)
p <- ggplot(Komen_top_hit, aes(x = Age, y = cg07303143)) + geom_point() 
p + geom_smooth(method='lm', se=TRUE) + theme_classic(base_size = 16) +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) +
  labs(list(x = "Subject Age", 
            y = paste(colnames(Komen_top_hit)[1], "(M-value)", sep=" ")))
dev.off()

# Examine distribution of P-values
hist(lf_Null$p.value[ , 'Age'])
hist(lf_Omega$p.value[ , 'Age'])
hist(lf_Null$p.value[ , 'BMI'])
hist(lf_Omega$p.value[ , 'BMI'])
hist(lf_Null$p.value[ , 'PregnantYes'])
hist(lf_Omega$p.value[ , 'PregnantYes'])

# Volcano plots for effect size
plot(lf_Omega$coefficients[ , 'Age'], -log10(lf_Omega$p.value[ , 'Age']), pch = 19)
plot(lf_Omega$coefficients[ , 'BMI'], -log10(lf_Omega$p.value[ , 'BMI']), pch = 19)
plot(lf_Omega$coefficients[ , 'PregnantYes'], -log10(lf_Omega$p.value[ , 'PregnantYes']), pch = 19)

# Adjustment for multiple comparisons (note: use of q-value with threshold of 0.01)
q_age_null <- qvalue(lf_Null$p.value[ , 'Age'], fdr.level = 0.01)
table(q_age_null$significant) # 787
q_age_null$pi0
qThresh = max(q_age_null$pvalues[q_age_null$qvalues<=0.01])
qThresh  # 0.000142765 
# Save output results from EWAS with age
Age_DMPs_null <- names(which(q_age_null$qvalues<=0.01)) # 4,099

# Adjustment for multiple comparisons (note: use of q-value with threshold of 0.01)
q_age <- qvalue(lf_Omega$p.value[ , 'Age'], fdr.level = 0.01)
table(q_age$significant) # 787
q_age$pi0
qThresh = max(q_age$pvalues[q_age$qvalues<=0.01])
qThresh  # 2.394087e-05 
# Save output results from EWAS with age
Age_DMPs <- names(which(q_age$qvalues<=0.01))
hyper_Age_DMPs <- names(which(lf_Omega$p.value[ , 'Age'] < 2.394087e-05 & lf_Omega$coef[ ,'Age'] > 0)) # 545
hypo_Age_DMPs <- names(which(lf_Omega$p.value[ , 'Age'] < 2.394087e-05  & lf_Omega$coef[ ,'Age'] < 0)) # 242

# Determine the CpGs that were no longer significant after adjustment.
Age_Unadjusted_Sig = Age_DMPs_null[-which(Age_DMPs%in%Age_DMPs_null)]
save(Age_Unadjusted_Sig, file="05.Genomic_Enrichment/Files/Komen_unsig_CpGs.RData")

# Assess risk factor relationship with DNAme genome-wide
q_bmi <- qvalue(lf_Omega$p.value[ , 'BMI'], fdr.level=0.01)
table(q_bmi$significant) # 0 significant CpGs
q_parity <- qvalue(lf_Omega$p.value[ , 'PregnantYes'], fdr.level = 0.01)
table(q_parity$significant) # 0 significant CpGs

# Store the total number of CpGs used in the Komen analysis
Komen_Analysis_CpGs <- rownames(lf_Omega$p.value) # All CpGs considered in analysis

# Vectors of CpG names will be used to determine enrichment/depletion in future analyses
save(list=c("Age_DMPs", "hyper_Age_DMPs", "hypo_Age_DMPs", "Komen_Analysis_CpGs"), file="05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")

##############################
# Create annotation for risk-factor related CpGs
##############################
library(GenomicRanges)
# Load annotation file that contains MAPINFO for CpG coordinates
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")

# Restrict the annotation file to only those CpGs used in the Komen analysis
annot_Komen_Age <- annot[Age_DMPs, ] 
annot_hyper_Age <- annot[hyper_Age_DMPs, ] 
annot_hypo_Age <- annot[hypo_Age_DMPs, ] 
annot_sub <- annot[Komen_Analysis_CpGs, ] # All CpGs considered in analysis

# Create start and stop genomic location indicators
annot_Komen_Age$hg19_start = as.numeric(annot_Komen_Age$MAPINFO)
# Create another new variable with MAPINFO for CpG end
annot_Komen_Age$hg19_end = as.numeric(annot_Komen_Age$MAPINFO)
#Create a new column that contans chromosome as 'chr' variable
annot_Komen_Age$chr = paste("chr", annot_Komen_Age$CHR, sep='')
# annot_sub_background <- annot_glioma[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg18_start', 'hg18_end', 'chr')]
annot_age_all <- annot_Komen_Age[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg19_start', 'hg19_end', 'chr')]
age_all_gr <- makeGRangesFromDataFrame(annot_age_all, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

# Create start and stop genomic location indicators
annot_hyper_Age$hg19_start = as.numeric(annot_hyper_Age$MAPINFO)
# Create another new variable with MAPINFO for CpG end
annot_hyper_Age$hg19_end = as.numeric(annot_hyper_Age$MAPINFO)
#Create a new column that contans chromosome as 'chr' variable
annot_hyper_Age$chr = paste("chr", annot_hyper_Age$CHR, sep='')
# annot_sub_background <- annot_glioma[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg18_start', 'hg18_end', 'chr')]
annot_age_hyper <- annot_hyper_Age[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg19_start', 'hg19_end', 'chr')]
age_hyper_gr <- makeGRangesFromDataFrame(annot_age_hyper, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

# Create start and stop genomic location indicators for hypomethylated with age
annot_hypo_Age$hg19_start = as.numeric(annot_hypo_Age$MAPINFO)
# Create another new variable with MAPINFO for CpG end
annot_hypo_Age$hg19_end = as.numeric(annot_hypo_Age$MAPINFO)
#Create a new column that contans chromosome as 'chr' variable
annot_hypo_Age$chr = paste("chr", annot_hypo_Age$CHR, sep='')
# annot_sub_background <- annot_glioma[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg18_start', 'hg18_end', 'chr')]
annot_age_hypo <- annot_hypo_Age[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg19_start', 'hg19_end', 'chr')]
age_hypo_gr <- makeGRangesFromDataFrame(annot_age_hypo, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

# 450K annotation file background
annot_sub$hg19_start = as.numeric(annot_sub$MAPINFO)
# Create another new variable with MAPINFO for CpG end
annot_sub$hg19_end = as.numeric(annot_sub$MAPINFO)
#Create a new column that contans chromosome as 'chr' variable
annot_sub$chr = paste("chr", annot_sub$CHR, sep='')
# annot_sub_background <- annot_glioma[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg18_start', 'hg18_end', 'chr')]
annot_sub_background <- annot_sub[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg19_start', 'hg19_end', 'chr')]
annot_background_gr <- makeGRangesFromDataFrame(annot_sub_background, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

#-------------------------------------------
# Convert GenomicRange objects to BED files
#-------------------------------------------
Age_All_df = data.frame(cbind(
  chrom = as.vector(seqnames(age_all_gr)),
  start = start(age_all_gr),
  end = end(age_all_gr)))
write.table(Age_All_df, file="05.Genomic_Enrichment/Files/Age_All_Komen.bed", quote=F, sep="\t", row.names=F, col.names=F)

Age_Hyper_df = data.frame(cbind(
  chrom=as.vector(seqnames(age_hyper_gr)),
  start=start(age_hyper_gr),
  end=end(age_hyper_gr)))
write.table(Age_Hyper_df, file="05.Genomic_Enrichment/Files/Age_Hyper_Komen.bed", quote=F, sep="\t", row.names=F, col.names=F)

Age_Hypo_df = data.frame(cbind(
  chrom=as.vector(seqnames(age_hypo_gr)),
  start=start(age_hypo_gr),
  end=end(age_hypo_gr)))
write.table(Age_Hypo_df, file="05.Genomic_Enrichment/Files/Age_Hypo_Komen.bed", quote=F, sep="\t", row.names=F, col.names=F)

# 5mC locations on the Illumina Chip
annot_back_df = data.frame(cbind(
  chrom=as.vector(seqnames(annot_background_gr)),
  start=start(annot_background_gr),
  end=end(annot_background_gr)))
write.table(annot_back_df, file="05.Genomic_Enrichment/Files/Annot_450k_Komen.bed", quote=F, sep="\t", row.names=F, col.names=F)

