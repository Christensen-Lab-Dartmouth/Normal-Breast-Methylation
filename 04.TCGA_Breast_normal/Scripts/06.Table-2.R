###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare Table for NDRI and TCGA Populations
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Include "count" function
library("plyr")

# Sara L defined function for covariate table building
Build_Table <- function(data_frame, variable_list){
  Table_DF <- rbind(c(paste("N =",nrow(data_frame)),"","",""),c("Variable","Mean (Range)","N (%)","NA's"))
  data_frame <- subset(data_frame,select=variable_list)
  for (col_num in 1:length(data_frame)){
    if (class(data_frame[,col_num]) != "factor"){
      mean_tmp <- round(mean(data_frame[,col_num],na.rm=T),digits=2)
      range_tmp <- round(range(data_frame[,col_num],na.rm=T),digits=2)
      if (any(is.na(data_frame[,col_num]))){
        NA_sum <- sum(is.na(data_frame[,col_num]))
      } else {
        NA_sum <- ""
      }
      Table_DF <- rbind(Table_DF,c(variable_list[col_num],
                                   paste(mean_tmp," (",range_tmp[1],"-",range_tmp[2],")",sep=""),
                                   "",NA_sum))
    } else {
      if (any(is.na(data_frame[,col_num]))){
        NA_sum <- sum(is.na(data_frame[,col_num]))
      } else {
        NA_sum <- ""
      }
      DF_tmp <- count(data_frame[,col_num])
      if (is.na(DF_tmp[nrow(DF_tmp),1])){
        DF_tmp <- DF_tmp[-nrow(DF_tmp),]
      }
      for (row in 1:nrow(DF_tmp)){
        if (row == 1){
          Table_DF <- rbind(Table_DF,c(variable_list[col_num],"","",NA_sum),
                            c(paste("  ",DF_tmp[row,1],sep=""),"",
                              paste(DF_tmp[row,2]," (",round((DF_tmp[row,2]/nrow(data_frame))*100,digits=2)," %)",sep=""),""))
        } else {
          Table_DF <- rbind(Table_DF,c(paste("  ",DF_tmp[row,1],sep=""),"",
                                       paste(DF_tmp[row,2]," (",round((DF_tmp[row,2]/nrow(data_frame))*100,digits=2)," %)",sep=""),""))
        }
      }
    }
    
  }
  return(Table_DF)
}

# Load Komen covariate file
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")
NDRI_covars <- read.csv("03.NDRI_Validation/Files/metaData-breast-160209-BSonly.csv", header = T, sep = ",", stringsAsFactors = F)
TCGA_covars <- read.csv("04.TCGA_Breast_Normal/Files/2016-08-04_TCGA_Normal_Breast_manifest.csv", header = T, sep = ",", stringsAsFactors = F)
TCGA_covars$Race[TCGA_covars$Race=='BLACK OR AFRICAN AMERICAN'] <- "AFRICAN AMERICAN"
TCGA_covars$Race[TCGA_covars$Race=='[Not Available]'] <- "NA"
TCGA_covars$Age = as.numeric(TCGA_covars$Age)
# Identify covariates of interest
NDRI_vars <- c("Sample_Age", "Sample_BMI")
TCGA_vars <- c("Age")

# Create Covariate Table
NDRI_T2 <- Build_Table(data_frame = NDRI_covars, variable_list = NDRI_vars)
TCGA_T2 <- Build_Table(data_frame = TCGA_covars, variable_list = TCGA_vars)

# Combine two populations for Table 2
validate_pops = rbind(NDRI_T2, TCGA_T2)
write.csv(validate_pops, "03.NDRI_Validation/Files/Validation-Populations.csv")


