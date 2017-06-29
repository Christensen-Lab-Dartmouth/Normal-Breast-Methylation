######################################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
######################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Enrichment using Fisher's Test for PCGTs and Enhancers from Roadmap to Epigenomics
# Only enhancer relationship is reported in the manuscript
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load annotation file that contains MAPINFO for CpG coordinates
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")
# PCGT data
load("05.Genomic_Enrichment/Files/PolycombComplete-120109.RData")
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")

# Columns need to be coded as.character to pull out desired information
annot$UCSC_RefGene_Accession <- as.character(annot$UCSC_RefGene_Accession)
annot$UCSC_RefGene_Name <- as.character(annot$UCSC_RefGene_Name)
annot$UCSC_RefGene_Group <- as.character(annot$UCSC_RefGene_Group)
tmp1 = strsplit(annot$UCSC_RefGene_Accession,";")
names(tmp1)=paste(1:length(tmp1),":",sep="")
tmp2 = unlist(tmp1)
tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
tmp4 = data.frame(UCSC_REFGENE_ACCESSION=tmp2,
                  rowid=tmp3, stringsAsFactors=FALSE)

#Need to curate Andy's PcG workspace before I proceed with this analysis of PcGs
tmp5 = merge(polycombTab,tmp4)
tmp6 = unique(tmp5$rowid)
annot$PcG = rep(0, dim(annot)[1])
annot$PcG[tmp6] = 1 
# save(annot, file="III.Risk_Factor_Associations/Files/Annot_PcG_1April2016.RData")

# Load processed annotation file to include 'PcG' column
load("05.Genomic_Enrichment/Files/Annot_PcG_1April2016.RData")
load("05.Genomic_Enrichment/Files/Komen_Analysis_CpGs.RData")

# Subset to CpGs considered in the analysis
annot_age = annot[Age_DMPs, ];table(annot_age$PcG)
annot_background = annot[Komen_Analysis_CpGs, ];table(annot_background$PcG)

Age_PCGT <- as.table(rbind(c(148, 49786),c(639, 340506)))
Fisher_PCGT <- fisher.test(Age_PCGT)
Fisher_PCGT  # 1.58 (CI:1.31 - 1.89)
Fisher_PCGT$p.value  # 1.740938e-06

##############################
# Roadmap to Epigenomics Data
##############################
# Step 1: Convert annotation file to GRanges object
library(GenomicRanges)
annot_background$hg19_start = as.numeric(annot_background$MAPINFO)
# Create another new variable with MAPINFO for CpG end
annot_background$hg19_end = as.numeric(annot_background$MAPINFO)
# Create a new column that contans chromosome as 'chr' variable
annot_background$chr = paste("chr", annot_background$CHR, sep='')
annot_background <- annot_background[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg19_start', 'hg19_end', 'chr')]
annot_background_gr <- makeGRangesFromDataFrame(annot_background, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")

# Step 2: Get the H3K4me1 (chromatin mark characterized by enhancer features)  ChIP-seq data
# from the Roadmap Epigenomics Project carried out on myoepithelial cells isolated from a 33 and 36 year old African American
# women's breast who was disease-free. Similar analyses could be completed with other
# available marks including: H3K4me3, H3K9ac, H3K9me2, H3K27me3, H3K36me3, and mRNAseq  

# Note on consolidated ChIP-seq data:
# Signal processing engine of MACSv2.0.10 peak caller to generate genome-wide signal coverage tracks
# Negative log10 of the Poisson p-value of ChIP-seq or DNase counts relative to expected background counts
# These signal confidence scores provides a measure of statistical significan of the observed enrichment.

# The sequencing data may come in BAM or .bigwig files. Need a package to load this data
library(rtracklayer)
sig.E027.H3K4me1.path = "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me1.pval.signal.bigwig"
# E027 = numeric epigenome identifier for consolidated data.
h3k4me1.E027 <- import.bw(sig.E027.H3K4me1.path)

#from the roadmap website: The -log10(palue) scores provide a convenient way to threshold signal (e.g. 2 corresponds to a a
# a p-value threshold of1e-2). A universal threshold of 2 provides good separation between signal and noise
sum(score(h3k4me1.E027)>2)
summary(score(h3k4me1.E027))
h3k4me1_sig_site = h3k4me1.E027[score(h3k4me1.E027)>=2, ]

median( width(h3k4me1_sig_site) )
#Plot the histogram
hist( width(h3k4me1_sig_site), nclass=25)
max( width(h3k4me1_sig_site) )

#We are interested in the overlap between the enhancer sites and 450K sites
# For the array
overlap_all_GR  <- subsetByOverlaps(annot_background_gr, h3k4me1_sig_site) 

#Determine which CpGs are in significant H3K4me1 sites
Breast_H3K4me1_CpGs = overlap_all_GR@elementMetadata$Name

# Create new column for Breast H3K4me1 sites (No=1) and (Yes=1)
annot_background$H3K4me1 = rep('0', dim(annot_background)[1])
annot_background$H3K4me1[rownames(annot_background)%in%Breast_H3K4me1_CpGs] = '1'

# Compare with conserved enhancers provided by Illumina annotation
table(annot_background$H3K4me1) # '1' = 89,731

# Reduce the background annotation to age-related DMPs
annot_age_enhancer <- annot_background[Age_DMPs, ]
table(annot_age_enhancer$H3K4me1) # '1' = 295

Age_H3K4 <- as.table(rbind(c(295, 89731),c(492, 300561)))
Fisher_H3K4 <- fisher.test(Age_H3K4)
Fisher_H3K4  # 2.00 (CI: 1.73 - 2.32)
Fisher_H3K4$p.value  # 7.138819e-20


