########################################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
######################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heatmap of age-related loci plotted alongside cell-type specific Roadmap to Epigenomics data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Roadmap to Epigenomics Data and Age-Related CpGs
load("05.Genomic_Enrichment/Files/Roadmap_Breast_Tissues.RData")
load("05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
sum(Age_DMPs%in%Roadmap_Breast_Tissues@rownames) # 691/787
Age_DMPs_Overlap <- Age_DMPs[which(Age_DMPs%in%Roadmap_Breast_Tissues@rownames)] # 691/787

Roadmap_dm <- data.matrix(cbind(as.numeric(Roadmap_Breast_Tissues@listData$Myoepithelial_Betas), 
                          as.numeric(Roadmap_Breast_Tissues@listData$Luminal_Betas), as.numeric(Roadmap_Breast_Tissues@listData$Adipose_Betas)))
colnames(Roadmap_dm) <- c("Myoepithelial", "Luminal", "Adipose")
rownames(Roadmap_dm) <- Roadmap_Breast_Tissues@rownames
# Only provide visualization for CpGs measured in both samples.
Roadmap_sub <- Roadmap_dm[Age_DMPs_Overlap, ]


# Load Betas that have been processed by minfi pipeline
Komen_Betas <- data.frame(fread("01.Preprocessing/Data/Komen_BMIQ_Betas_17May2016.csv"), row.names=1)
Betas = Komen_Betas[ ,order(colnames(Komen_Betas), decreasing=T)]
colnames(Betas) <- substring(colnames(Betas), 2, nchar(colnames(Betas))[1])

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

Betas_Komen_All_Age <- Betas_Komen[Age_DMPs, ]
Betas_Komen_sub <- Betas_Komen[Age_DMPs_Overlap, ]

Komen_Ref_CellTypes <- merge(Betas_Komen_sub, Roadmap_sub, by="row.names")
rownames(Komen_Ref_CellTypes) <- Komen_Ref_CellTypes[,1]
Komen_Ref_CellTypes <- Komen_Ref_CellTypes[,-1]
########################
# Load annotation file that contains MAPINFO for CpG coordinates
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")

# Restrict the annotation file to only those CpGs used in the Komen analysis
annot_Road <- annot[rownames(Komen_Ref_CellTypes), ] 
levels(annot_Road$Relation_to_UCSC_CpG_Island)[levels(annot_Road$Relation_to_UCSC_CpG_Island)==""] <- "Open Sea"

# Collapse "North" and "South" nomenclature for CpG islands
CpGs_Annotated <- gsub("^[N|S]_","",annot_Road$Relation_to_UCSC_CpG_Island) 

# Column colors for  CpG island regions
colregions <- c("black", "blue", "ghostwhite", "darkgray")
# Get CpG Island Colors
Regions <- c("1","2","3","4")
colorRegionsVector <- c(as.character(as.numeric(as.factor(CpGs_Annotated))))
for(i in 1:length(colregions))
{
  colorRegionsVector[colorRegionsVector == Regions[i]] <- colregions[i]
}

row_colors = as.matrix(t(colorRegionsVector))

#################
###Row colors ###
#################
# Ever pregnant?
PREGNANT <- c('Yes', 'No')
Pregcolors <- c("green", "white")
#Get Subtype Colors
colorPREGVector <- covars_Komen$Pregnant
for(i in 1:length(Pregcolors))
{
  colorPREGVector[colorPREGVector == PREGNANT[i]] <- Pregcolors[i]
}

#Age at donation  
colorAge <- hsv(0.01, covars_Komen$Age/max(covars$Age), 1)

#BMI at donation 
colorBMI <- hsv(0.01, covars_Komen$BMI/max(covars$BMI), 1)

# Pink = myoepithelial, blue = luminal, adipose = yellow
Roadmap_colors <- c("pink","blue", "yellow")

# columncolors <- as.matrix(t(rbind(c(colorPREGVector,Roadmap_colors), c(colorAge, Roadmap_colors), c(colorBMI, Roadmap_colors))))
columncolors <- as.matrix(t(rbind(c(colorAge, Roadmap_colors))))
colnames(columncolors) = "Subject Age"
col_colors = as.matrix(columncolors[1:100, ])
# colnames(columncolors) <- c("Ever Pregnant", "Subject BMI", "Subject Age")

heatcols <- colorRampPalette(c("yellow","blue"))(n = 600)
colbreaks <- seq(0,1, 1/600)

####################
source("/Users/kevinjohnson/Komen_Normal_5mC/02.RefFreeEWAS/Scripts/07.make_heatmaps.R")
####################
png("05.Genomic_Enrichment/Figures/Komen_Roadmap_Heatmap.png", height=7*300, width=7*300, res=300)
heatmap.3(as.matrix(Komen_Ref_CellTypes[ , 1:100]), col = heatcols, breaks = colbreaks, trace = "none", labCol = F, cexCol = 0.5, key = F, density.info = "none", 
          scale = "none", dendrogram='column', main = "", ColSideColors = col_colors, RowSideColors = row_colors, ColSideColorsSize = 2,
          cexRow = 1.5, labRow = F)
dev.off()

png("05.Genomic_Enrichment/Figures/Komen_Roadmap_CpG_Legend.png", height=7*300, width=7*300, res=300)
plot(0,type='n',axes=FALSE,ann=FALSE)
par(lend = 1) 
legend("top",      
       legend = c("CpG Island", "Shore", "Shelf", "Open Sea"),
       col = c("black", "darkgray", "ghostwhite",  "blue"), 
       bty ='n',
       lwd = 5)
dev.off()

white_mat <- matrix(0.5, nrow = 691, ncol = 100)
HM_Komen <- heatmap.3(as.matrix(Komen_Ref_CellTypes[ , 1:100]))
Road_Colors <- as.matrix(c(rep("white",100), "pink", "blue", "yellow"))
Roadmap_colors <- c("pink","blue", "yellow")

png("05.Genomic_Enrichment/Figures/Roadmap_Clustering.png", height=7*300, width=7*300, res=300)
heatmap.3(cbind(white_mat[ ,], Komen_Ref_CellTypes[rev(HM_Komen$rowInd), 101:103]), dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none',
          labRow = F, col = heatcols, key = F, ColSideColors = Road_Colors, ColSideColorsSize = 2, labCol=F)
dev.off()

