###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Alex Titus, Jess King, Andres Houseman, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Use the Whole Genome Bisulfite Sequencing Data from Roadmap to convert to 450k
# Use approach presented in Titus et al Bioninformatics 2016
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Location of WGBS data files processed by Alex Titus
setwd("/Users/kevinjohnson/Documents/WGBS_450K/Roadmap_samples/")

# From the Roadmap to Epigenomics (GEO)
load("WGBSliftover_Roadmap_BRCA_GSM1127059_myo.rdata")
myoepithelial <- TCGA_WGBSto450K
myoepithelial = myoepithelial[ ,c(5,9)]
colnames(myoepithelial) = c("Myoepithelial_Betas", "Island_Status") 
rm(TCGA_WGBSto450K)

load("WGBSliftover_Roadmap_BRCA_GSM1127125_lum.rdata")
luminal <- TCGA_WGBSto450K
luminal = luminal[ ,c(5,9)]
colnames(luminal) = c("Luminal_Betas", "Island_Status") 
rm(TCGA_WGBSto450K)

load("WGBSliftover_Roadmap_ADIP_GSM1120331.rdata")
adipose <- TCGA_WGBSto450K
adipose = adipose[ ,c(5,9)]
colnames(adipose) = c("Adipose_Betas", "Island_Status") 
rm(TCGA_WGBSto450K)

# Common set of CpGs 
breast_cells = merge(myoepithelial, luminal, by="row.names")
rownames(breast_cells) = breast_cells$Row.names; breast_cells$Row.names = NULL
major_breast_cells = merge(breast_cells, adipose, by="row.names")
rownames(major_breast_cells) = major_breast_cells$Row.names; major_breast_cells$Row.names = NULL
Roadmap_Breast_Tissues = major_breast_cells[ , c("Myoepithelial_Betas", "Luminal_Betas", "Adipose_Betas")]
save(Roadmap_Breast_Tissues, file="Roadmap_Breast_Tissues.RData") # 409,297 CpGs 

# Fibroblast
load("fibroRRBSliftover450k.RData")
fibroblast <- WGBSliftover
fibroblast <- as.data.frame(fibroblast)
rownames(fibroblast) = fibroblast$probe
rm(WGBSliftover)
fibroblast = fibroblast[ ,c("V5", "Relation_to_Island")]
colnames(fibroblast) = c("Fibroblast_Betas", "Island_Status") 
Roadmap_Breast_4Cells = merge(Roadmap_Breast_Tissues, fibroblast, by="row.names")
rownames(Roadmap_Breast_4Cells) = Roadmap_Breast_4Cells $Row.names; Roadmap_Breast_4Cells$Row.names = NULL
save(Roadmap_Breast_4Cells, file="Roadmap_Breast_4Cells.RData") # 90, 668 CpGs 


