#######################################
#  Generate Heatmap for age-related   #
#######################################
library(data.table)
library(gplots)
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

# All age-related DMPs from Komen
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")

# Subet to only those age-related CpGs
Komen_Age_Related = Betas_Komen[Age_DMPs, ]

# Use custom heatmap.3 function
source("/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Scripts/07.make_heatmaps.R")

# Load annotation file that contains MAPINFO for CpG coordinates
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")

# Restrict the annotation file to only those CpGs used in the Komen analysis
annot_Komen_Age <- annot[Age_DMPs, ] 

library(plyr)
revalue(annot_Komen_Age$Relation_to_UCSC_CpG_Island, c("Open Sea"="Open Sea", "N_Shore"="Shore", "S_Shore"="Shore",
                                                       "S_Shelf"= "Shelf", "N_Shelf"="Shelf", "Island"="Island"))
levels(annot_Komen_Age$Relation_to_UCSC_CpG_Island)[levels(annot_Komen_Age$Relation_to_UCSC_CpG_Island)==""] <- "Open Sea"
# Collapse "North" and "South" nomenclature for CpG islands
CpGs_Annotated <- gsub("^[N|S]_","",annot_Komen_Age$Relation_to_UCSC_CpG_Island) 

# Age at donation  
colorAge <- hsv(0.01, covars_Komen$Age/max(covars_Komen$Age), 1)

# BMI at donation 
colorBMI <- hsv(0.1, covars_Komen$BMI/max(covars_Komen$BMI), 1)

# Column colors for  CpG island regions
colregions <- c("black", "blue", "ghostwhite", "darkgray")
# Get CpG Island Colors
Regions <- c("1","2","3","4")
colorRegionsVector <- c(as.character(as.numeric(as.factor(CpGs_Annotated))))
for(i in 1:length(colregions))
{
  colorRegionsVector[colorRegionsVector == Regions[i]] <- colregions[i]
}

#Define the color palette
heatcols <- colorRampPalette(c("yellow","black",'blue'))(n = 1000)
png("/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Figures/Age_Related_CpGs_Heatmap.png", height=7*300, width=7*300, res=300)
heatmap.2(Komen_Age_Related, col=heatcols, key.title="", key.xlab="DNAme", dendrogram='none', 
          RowSideColors =  colorRegionsVector, ColSideColors=colorAge, main = "", margins=c(5,8), 
          xlab="", labRow=FALSE, labCol=FALSE, trace = "none")

par(lend = 1) 
legend("top",      
       legend = c("CpG Island", "Shore", "Shelf", "Open Sea"),
       col = c("black", "darkgray", "ghostwhite",  "blue"), 
       bty ='n',
       lwd = 5)
dev.off()



