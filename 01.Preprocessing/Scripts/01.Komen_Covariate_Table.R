###########################
#  DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
# Adapted from: Sara Lundgren (30 March 2015)
######################################################
# ~~~~~~~~~~~~~~~~~~
# Generate a covariate table from the Komen Tissue Bank subjects
# ~~~~~~~~~~~~~~~~~~
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
covariates <- read.csv("01.Preprocessing/Files/Komen_Manifest_20July2015.csv", header = T, sep = ",")

# Identify covariates of interest
Komen_vars <- c("Age", "BMI", "Pregnant", "Menstrual.Status", "BLDCAN", "Race", "Gail_Score")

# Create Covariate Table
Komen_T1 <- Build_Table(data_frame = covariates, variable_list = Komen_vars)

# Write table out
write.csv(Komen_T1, "01.Preprocessing/Files/Komen_Table1.csv")
