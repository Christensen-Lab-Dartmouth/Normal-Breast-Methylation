###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Test the relation between TCGA normal breast DNAme and gene expression for age-related DMPs
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# NECESSARY libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# All age-related DMPs from Komen
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")

# Load Betas that have been processed by minfi
TCGA_Betas <- data.frame(fread("04.TCGA_Breast_Normal/Files/TCGA_BMIQ_Betas_25July2016.csv"), row.names=1)
colnames(TCGA_Betas) <- substring(colnames(TCGA_Betas), 2, length(colnames(TCGA_Betas)))
Betas = TCGA_Betas[ ,order(colnames(TCGA_Betas), decreasing=T)]

# Load covariate file
covariates <- read.csv("04.TCGA_Breast_Normal/Files/2016-08-04_TCGA_Normal_Breast_manifest.csv", header = T, sep = ",", stringsAsFactors = F)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars_TCGA = covariates[order(covariates$Match_IDs, decreasing=T),]

# Use only those CpGs that were associated with age in normal breast tissue
Betas_TCGA <- data.matrix(Betas[Age_DMPs, ])
all(names(Betas_TCGA) == covars_TCGA$Match_IDs)
# Provide new sample TCGA sample identifiers
colnames(Betas_TCGA) <- covars_TCGA$Sample_Name

# Load TCGA normal expression data
load("07.Gene_Expression/Files/TCGA_BRCA_Normal_RNAseq_rsem.RData")
rownames(NormBRCARNAseq) <- substring(rownames(NormBRCARNAseq), 1, 12)

# Remove additional information from sample name in the expression data set to match with DNAme
genes <- colnames(NormBRCARNAseq)
tmp <- sapply(genes, function(x) unlist(strsplit(x, "[|]"))) 
tmp2 <- as.character(tmp[1, ])
colnames(NormBRCARNAseq) <- tmp2

# Use only those samples that have both DNAme and Gene expression in TCGA
Expression_match <- colnames(Betas_TCGA)[colnames(Betas_TCGA)%in%rownames(NormBRCARNAseq)]
Betas_TCGA_matched <- Betas_TCGA[ , Expression_match]
RNAseq_TCGA_matched <- NormBRCARNAseq[Expression_match, ]
RNAseq_TCGA_matched <- t(RNAseq_TCGA_matched)
all(colnames(Betas_TCGA_matched)==colnames(RNAseq_TCGA_matched))


# Perform correlation analyses using the gene name. CpGs without gene name get dropped.
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")
annot$UCSC_RefGene_Name <- as.character(annot$UCSC_RefGene_Name)
annot$UCSC_RefGene_Group <- as.character(annot$UCSC_RefGene_Group)
annotSel <- annot[Age_DMPs, ]

# The process flow for the loop is:
# 1. Pull out CpG name 
# 2. Find the gene associated with this CpG 3. 
# 3. If the gene name is NA or blank. Going to skip those
# 4. Find out where in the RNA Seq file the Gene in located (there are no duplicate gene names)
# 5. If there is nothing there, skip it
# 6. Calculate the Pearson and Spearman correlations on untransformed values
# 7. past each line in the .txt file over and over... and over (etc.) until the last line is read

# Transpose the betas from the methylation object so that they have the same orientation as the RNASeq file. Double check to make sure that the sample IDs have the same orientation
cpg = Betas_TCGA_matched
CpGgenes <- annotSel
CpGgenes$UCSC_RefGene_Name = as.character(CpGgenes$UCSC_RefGene_Name)

#Listing the output file where R will dump the results
myoutf1 = "/Users/kevinjohnson/Komen_Normal_5mC/07.Gene_Expression/Files/Normal_RNAseq_TCGA_BRCA_ageDMPs.txt"

conOut = file(myoutf1, "w")

for(h in 1:nrow(cpg)){
  mycpg = row.names(cpg)[h]
  xx = CpGgenes[mycpg, 21]
  xx = unlist(strsplit(xx, ";"))[1]
  if(is.na(xx)) next  	
  if(xx=="") next	
  se = which(rownames(RNAseq_TCGA_matched)==xx)
  if(length(se)==0) next
  corr = cor(as.numeric(RNAseq_TCGA_matched[se,]), as.numeric(cpg[h,]), method="s")	
  if(is.na(corr)) next
  corr2 = cor.test(as.numeric(RNAseq_TCGA_matched[se,]), as.numeric(cpg[h,]), method="s")	
  curLine = paste(mycpg, CpGgenes[mycpg,21], CpGgenes[mycpg,23], round(corr2$estimate, 6), corr2$p.value, sep="\t")
  writeLines(curLine, conOut)	
}
close(conOut)


# Generate figures to support relationship with expression
Meth.Express <- read.table("/Users/kevinjohnson/Komen_Normal_5mC/07.Gene_Expression/Files/Normal_RNAseq_TCGA_BRCA_ageDMPs.txt", sep = "", stringsAsFactors = F)
colnames(Meth.Express) <- c("CpG_ID", "Gene_Name", "Gene_Region", "Correlation_Coef", "P_value")
png("/Users/kevinjohnson/Komen_Normal_5mC/07.Gene_Expression/Figures/5mC-Express-Breast-Normal-Histogram.png", width=7*300, height=7*300, res=300)
hist(Meth.Express$P_value, main="5mC-Expression Correlations TCGA Normal Breast\n Age DMPs Komen", xlab="P-value")
dev.off()

# How many CpGs, that track to genes, are associated with gene expression at nominal level
sum(Meth.Express$P_value<0.05) # 259/630

# Generate Gene Region Plot of 5mC and gene expresison from 84 normal breast tissues in the TCGA
Meth.Express$iso_gene_names = sapply(strsplit(Meth.Express$Gene_Name, ";"), '[[', 1)
Meth.Express$iso_gene_regions = sapply(strsplit(Meth.Express$Gene_Region, ";"), '[[', 1)
# Factor to relevel
Meth.Express$iso_gene_regions <- factor(Meth.Express$iso_gene_regions, levels = c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body","3'UTR"))
Meth.Express.sig <- Meth.Express[Meth.Express$P_value<0.05, ]
Meth.Express.sig$P_value[Meth.Express.sig$P_value==0] = 2.2e-16
# Make sure the bubbles don't sit on top of one another
jit <- position_jitter(width = .4)
# Plot for the correlation, significance, and region of DNA methylation/expression changes
# Layer graphics for each gene region
p <- ggplot(Meth.Express.sig, aes(x=iso_gene_regions, y=Correlation_Coef, size= -log10(P_value), color=iso_gene_regions, alpha=1/5)) + ylim(-1,1) + geom_jitter(position = jit)

# Plot the gene names for the 25 CpGs with the lowest P-value
png("/Users/kevinjohnson/Komen_Normal_5mC/07.Gene_Expression/Figures/5mC-Express-Breast-Normal-SigOnly-Figure.png", width=7*300, height=7*300, res=300)
p + xlab("Genomic Regions") +ylab("Correlation Coefficient")  +  geom_text_repel(data=filter(Meth.Express.sig, P_value<4.550339e-08), aes(label=iso_gene_names)) +
  theme_bw() + ggtitle("5mC-Expression Correlations") + 
  geom_hline(yintercept = 0) + scale_alpha(guide='none') 
dev.off()

