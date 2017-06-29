########################################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
######################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Locus Overlap Analysis of age-related CpGs 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code below is based on instructions written by: Nathan Sheffield
# source("https://bioconductor.org/biocLite.R"); biocLite("LOLA")
# See R citation at: browseVignettes("LOLA")
# See Fig.1 of Sheffield & Bock (2015) for examples of ranked LOLA enrichment results
# Project directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
# Necessary libraries
library(LOLA)
library(GenomicRanges) # LOLA requires loading of datasets as GRanges object

# Define the universe & the gene set (the concept of "universe" is well-explained in Sheffield's "Choosing a LOLA Universe")
array <- readBed('05.Genomic_Enrichment/Files/Annot_450k_Komen.bed') # File for universe 
age_related_all <- readBed('05.Genomic_Enrichment/Files/Age_All_Komen.bed') # Data for query input
age_related_hyper <- readBed('05.Genomic_Enrichment/Files/Age_Hyper_Komen.bed') # Data for query input
age_related_hypo <- readBed('05.Genomic_Enrichment/Files/Age_Hypo_Komen.bed') # Data for query input

# Set Universe as CpGs considered in Komen analysis
universe <- GRanges(array)
geneset_hyper <- GRanges(age_related_hyper)
geneset_hypo <- GRanges(age_related_hypo)

# Compute gene enrichment set (see document "Using LOLA Core")
# Load the LOLA core (cached version; hg19/38, ENCODE TFBS, UCSC CGIs, Citrome epigenome)
library(devtools)
install_github("sheffien/simpleCache")

lolaDB <- loadRegionDB("/Users/kevinjohnson/Komen_Normal_5mC/05.Genomic_Enrichment/LOLACore/hg19/")
locResults_hyper <- runLOLA(geneset_hyper, universe, lolaDB)
dim(locResults_hyper)
locResults_hyper[1:10,] #list the top 10 results
locResults_hyper[order(meanRnk, decreasing=FALSE),][1:20,]
locResults_hyper[order(maxRnk, decreasing=FALSE),][1:20,]
locResults_hyper[collection=="sheffield_dnase", ][order(meanRnk, decreasing=FALSE),]
locResults_hyper[collection=="encode_tfbs", ][order(meanRnk, decreasing=FALSE),]
locResults_hyper[collection=="codex", ][order(meanRnk, decreasing=FALSE),]
locResults_hyper[collection=="cistrome_epigenome", ][order(meanRnk, decreasing=FALSE),]
locResults_hyper[collection=="ucsc_features", ][order(meanRnk, decreasing=FALSE),]

# There are 61 CpGs that lose DNA methylation with age. Test for enrichment separately.
locResults_hypo <- runLOLA(geneset_hypo, universe, lolaDB)
locResults_hypo[1:10,] #list the top 10 results
locResults_hypo[collection=="encode_tfbs", ][order(meanRnk, decreasing=FALSE),]

# Output of all enrichment files
writeCombinedEnrichment(locResults_hyper, outFolder= "05.Genomic_Enrichment/lolaResults_hyper", includeSplits=TRUE)
writeCombinedEnrichment(locResults_hypo, outFolder= "05.Genomic_Enrichment/lolaResults_hypo", includeSplits=TRUE)

