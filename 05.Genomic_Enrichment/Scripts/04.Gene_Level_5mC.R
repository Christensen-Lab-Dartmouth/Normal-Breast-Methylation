###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Provide CpG gene list by region
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
load("/Users/kevinjohnson/Komen_Normal_5mC/01.Preprocessing/Data/annot.RData")

# Restrict the annotation file to the age-related CpGs in Komen
load("/Users/kevinjohnson/Komen_Normal_5mC/05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
annotSel <- annot[Age_DMPs, ] 

#Forgot to set stringsAsFactors = F in previous 'annot' data.frame
annotSel$UCSC_RefGene_Name <- as.character(annotSel$UCSC_RefGene_Name)
annotSel$UCSC_RefGene_Group <- as.character(annotSel$UCSC_RefGene_Group)

# Get an index of CpGs by gene region
tmp1 = strsplit(annotSel$UCSC_RefGene_Name,";")
names(tmp1)=paste(1:length(tmp1),":",sep="")
tmp2 = unlist(tmp1)
tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
tmp4 = strsplit(annotSel$UCSC_RefGene_Group,";")
names(tmp4)=paste(1:length(tmp4),":",sep="")
tmp5 = unlist(tmp4)
tmp6 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp5)))
all(tmp3==tmp6)

GeneAnnotation = unique(data.frame(UCSC_REFGENE_NAME=tmp2,UCSC_REFGENE_GROUP=tmp5,
                                   rowid=tmp3, stringsAsFactors=FALSE))
GeneIndex = split(GeneAnnotation$rowid,
                  with(GeneAnnotation,paste(UCSC_REFGENE_NAME,UCSC_REFGENE_GROUP,sep=":")))
GeneIndexN = sapply(GeneIndex, length)

tmp1 = sapply(GeneIndex, function(u) u[which.min(annotSel$MAPINFO[u])])
tmp2 = sapply(GeneIndex, function(u) u[which.max(annotSel$MAPINFO[u])])
GeneSym1 = annotSel$UCSC_RefGene_Name[tmp1]
GeneSym2 = annotSel$UCSC_RefGene_Name[tmp2]

tmp = strsplit(names(GeneIndexN),":")
GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),
                         Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN,
                         stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
write.csv(GeneResults, file="05.Genomic_Enrichment/Files/GeneSummary-Age-Related-Komen.csv", row.names=FALSE)

#######################
# Background of annotation file
#######################
load("05.Genomic_Enrichment/Files/Komen_Analysis_CpGs.RData")
annotSel <- annot[Komen_Analysis_CpGs, ] 

#Forgot to set stringsAsFactors = F in previous 'annot' data.frame
annotSel$UCSC_RefGene_Name <- as.character(annotSel$UCSC_RefGene_Name)
annotSel$UCSC_RefGene_Group <- as.character(annotSel$UCSC_RefGene_Group)

# Get an index of CpGs by gene region
tmp1 = strsplit(annotSel$UCSC_RefGene_Name,";")
names(tmp1)=paste(1:length(tmp1),":",sep="")
tmp2 = unlist(tmp1)
tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
tmp4 = strsplit(annotSel$UCSC_RefGene_Group,";")
names(tmp4)=paste(1:length(tmp4),":",sep="")
tmp5 = unlist(tmp4)
tmp6 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp5)))
all(tmp3==tmp6)

GeneAnnotation = unique(data.frame(UCSC_REFGENE_NAME=tmp2,UCSC_REFGENE_GROUP=tmp5,
                                   rowid=tmp3, stringsAsFactors=FALSE))
GeneIndex = split(GeneAnnotation$rowid,
                  with(GeneAnnotation,paste(UCSC_REFGENE_NAME,UCSC_REFGENE_GROUP,sep=":")))
GeneIndexN = sapply(GeneIndex, length)

tmp1 = sapply(GeneIndex, function(u) u[which.min(annotSel$MAPINFO[u])])
tmp2 = sapply(GeneIndex, function(u) u[which.max(annotSel$MAPINFO[u])])
GeneSym1 = annotSel$UCSC_RefGene_Name[tmp1]
GeneSym2 = annotSel$UCSC_RefGene_Name[tmp2]

tmp = strsplit(names(GeneIndexN),":")
GeneResults_background = data.frame(Gene=sapply(tmp,function(u)u[1]),
                         Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN,
                         stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)

Background_Check <- merge(GeneResults, GeneResults_background, by="row.names")
write.csv(Background_Check, file="05.Genomic_Enrichment/Files/GeneSummary-Age-Related-Komen-Background-Check.csv", row.names=FALSE)
