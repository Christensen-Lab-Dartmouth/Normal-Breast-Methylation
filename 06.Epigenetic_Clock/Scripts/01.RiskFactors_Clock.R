###########################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
###########################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test the association of breast cancer risk factors and age acceleration
# as calculated by the epigenetic clock (http://labs.genetics.ucla.edu/horvath/dnamage/)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########################################
# Part I: Data cleaning
#########################################
library(ggplot2)
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

# Load Epigenetic Clock output and determine relation between DNAme-Chronological Age
covariates <- read.csv("01.Preprocessing/Files/Komen_Manifest_20July2015.csv", header = T, sep = ",", stringsAsFactors = F,row.names=1)
covariates$Match_IDs <- paste(covariates$Sentrix_ID, covariates$Sentrix_Position, sep="_")
covars = covariates[order(covariates$Match_IDs, decreasing=T),]
full_covariates <- read.csv("01.Preprocessing/Files/2016_08_08_Full_Covariate_Komen_Samples.csv", header=T, row.names=1)
tmp_covar <- merge(covars, full_covariates, by="row.names")
rownames(tmp_covar) <- tmp_covar$Row.names
tmp_covar$Row.names <- NULL
tmp_covar$Komen_ID <- rownames(tmp_covar)

# Results from submitted 450K array data (unadjusted for BMIQ as Horvath implements normalizatio procedure)
Clock_output <- read.csv("06.Epigenetic_Clock/Files/Betas_Clock_K.output.csv", header=TRUE)
Horvath_clock_plot <- ggplot(Clock_output, aes(x = Age, y = DNAmAge)) + geom_point() 
Horvath_clock_plot + geom_smooth(method='lm', se=TRUE) 

# Properly orient the data
covar_clock <- merge(tmp_covar, Clock_output, by.x="Match_IDs", by.y="X")
rownames(covar_clock) <- covar_clock$Row.names
covar_clock$Row.names <- NULL

################################################
#Part 2: DNAme-Chronological Age Relationships
################################################
# Plot the residuals of DNA methylation age regressed on chronological age
Clock_regress <- lm(DNAmAge~Age, data=Clock_output)
summary(Clock_regress)$coefficients[2,4] # Exact P-value, 4.08E-51
Clock_resid <- resid(Clock_regress)
png("06.Epigenetic_Clock/Figures/Residuals_DNAme_Age_Komen.png", height=7*300, width=7*300, res=300)
plot(Clock_output$Age, Clock_resid, ylab="Residuals", xlab="Chronological Age", main="Komen Normal (n=100)", pch=19)
abline(0, 0)
dev.off()

# The AgeAccel variable is normally distributed
hist(covar_clock$AgeAccelerationResidual, main="Histogram of Age Acceleration", xlab='Age Acceleration Residual')

# Produce a ggtheme for the background of the plot
covar_clock$Race.x[covar_clock$Race.x=="AFRNAMER"] = "AFRICAN-AMERICAN"
Breast_Clock <- ggplot(covar_clock, aes(x=Age.x, y= DNAmAge, color= Race.x,  size=abs(AgeAccelerationResidual) ) )  + geom_point(alpha=1/3) +  ylim(10,100) + xlim(10,100) + 
  geom_abline(intercept= 0, slope=1, colour='black', alpha=1/3) + scale_size(range = c(2.5,7.5)) + theme_bw() 
png("06.Epigenetic_Clock/Figures/AgeAccel-Normal-Figure.png", height=7*300, width=7*300, res=300)
Breast_Clock + annotate('text',x=70, 15, label='DNAme-Age P = 4.1E-51', alpha=0.8) + xlab("Chronological Age")  + ylab("DNA Methylation Age") +
  labs(size="Age-Acceleration (Years)", color="Race") + theme(legend.key = element_blank())
dev.off()

############################################
# Part 3: Age-Acceleration and Risk Factors
############################################
# Univariate analyses for each of the major risk factors in the Komen data set
# Create BMI bins for overweight (BMI greater than 25) and normal/underweight (< 25)
covar_clock$BMI_group <- ifelse(as.numeric(covar_clock$BMI.x>=25.0), "overweight", "normal") 
table(covar_clock$BMI_group)
# Treat BMI as a continuous and categorical variable
BMI_num_accel <- glm(as.numeric(AgeAccelerationResidual)~as.numeric(BMI.x),  family= gaussian, data = covar_clock)
BMI_cat_accel <- glm(as.numeric(AgeAccelerationResidual)~as.factor(BMI_group),  family= gaussian, data = covar_clock)
summary(BMI_num_accel)
summary(BMI_cat_accel)

# Gail Score for those women over age 35
BRCAT_accel <- glm(as.numeric(AgeAccelerationResidual)~as.numeric(Gail.Score),  family= gaussian, data = covar_clock)
summary(BRCAT_accel) # P = 0.62

# Age.at.Menarche
Men_accel <- glm(as.numeric(AgeAccelerationResidual)~as.numeric(Age.at.Menarche),  family= gaussian, data = covar_clock)
summary(Men_accel) # P = 0.56

# Change naming scheme for drinks per week
covar_clock$Drinks.per.Week <- as.factor(covar_clock$Drinks.per.Week)
covar_clock$drinking_status <- ifelse(covar_clock$Drinks.per.Week == " ", 0,  ifelse(covar_clock$Drinks.per.Week == " < 7", 1, 2) )
Drinks_per_week <- glm(as.numeric(AgeAccelerationResidual) ~ as.factor(drinking_status),  family= gaussian, data = covar_clock)
summary(Drinks_per_week) # > 7 drinks per week. P = 0.39

# Parity/Pregnancy 
Pregnancy <- glm(as.numeric(AgeAccelerationResidual)~as.factor(Pregnant), family= gaussian, data = covar_clock)
summary(Pregnancy)# P = 0.99

# Family History (A Blood relative diagnosed with breast or ovarian cancer)
FamHx <- glm(as.numeric(AgeAccelerationDiff)~ as.factor(BLDCAN), family= gaussian, data = covar_clock)
summary(FamHx) # P = 0.37

# Race (White, Hispanic, and African-American)
covar_clock$Race.x <- relevel(as.factor(covar_clock$Race.x), ref = "WHITE") # n = 86 white women
Race <- glm(as.numeric(AgeAccelerationResidual)~ as.factor(Race.x), family = gaussian, data = covar_clock)
summary(Race) # African-American P = 0.035*, Hispanic P = 0.81

# Age at First Birth
Age_at_first_birth <- glm(as.numeric(AgeAccelerationResidual)~Age.at.First.Birth,  family= gaussian, data = covar_clock)
summary(Age_at_first_birth) # P = 0.32

# Height
covar_clock$total_height <- (covar_clock$Feet*12+covar_clock$Inches)
hist(covar_clock$total_height)
Female_Height <- glm(as.numeric(AgeAccelerationResidual)~as.numeric(total_height),  family= gaussian, data = covar_clock)
summary(Female_Height) # P = 0.42

##### Multivariate regression #####
# DNAm age was regressed on age, BMI, and other covariates
# Lifestyle factors: BMI, Children, Drinks per week
Multi_var_lifestyle <- glm(as.numeric(AgeAccelerationResidual)~as.factor(Ever.Smoked.)+as.numeric(BMI.x)+as.factor(Pregnant)+as.factor(drinking_status), 
                 data = covar_clock, family = gaussian)
summary(Multi_var_lifestyle)

# Biological factors: Height, Family History, Age at Menarche, and Race
Multi_var_biological <- glm(as.numeric(AgeAccelerationResidual)~as.factor(BLDCAN)+as.numeric(Age.at.Menarche)+as.factor(Race.x), 
                 data = covar_clock, family = gaussian)
summary(Multi_var_biological)
table(covar_clock$Race.x)

# Complete Model: BMI*Menopause interaction, Parity, Family History, Age at Menarche, and Race
covar_clock$Menstrual.Status.x <- relevel(as.factor(covar_clock$Menstrual.Status.x), ref = "Pre-Menopausal")
Multi_var_full_model <- lm(as.numeric(AgeAccelerationResidual)~as.factor(BMI_group)+as.factor(BLDCAN)+as.factor(Pregnant)+as.factor(Race.x)+as.factor(drinking_status), 
                            data = covar_clock)
summary(Multi_var_full_model) # B = 4.0, P = 0.035 for African-American AgeAccel. 
summary(Multi_var_full_model_interact) # B = 3.8785, P = 0.0292 for African-American AgeAccel.
interact_fit <- anova(Multi_var_full_model, Multi_var_full_model_interact)
summary(interact_fit)

# Examine the distribution
aggregate(AgeAccelerationResidual~BMI_group+Menstrual.Status.x,data = covar_clock, median)
table(covar_clock$BLDCAN)

# Restrict the covariates to WHITE women only (the majority)
covar_clock_WHITE <- covar_clock[covar_clock$Race.x == "WHITE", ]
Multi_var_WHITE <- glm(as.numeric(AgeAccelerationResidual)~as.numeric(BMI.x)*as.factor(Menstrual.Status.x)+as.factor(Pregnant)+as.factor(BLDCAN)+as.numeric(Age.at.Menarche), 
                 data= covar_clock_WHITE, family= gaussian)
summary(Multi_var_WHITE)

############################################################
# Part 4: Quantitate the overlap between Horvath and age-related loci
############################################################
Horvath_clock_annot <- read.csv("/Users/kevinjohnson/Komen_Normal_5mC/06.Epigenetic_Clock/Files/GB_Horvath.csv", header=T, row.names=1)
Horvath_CpGs = rownames(Horvath_clock_annot)
# Inspect how many age-related breast CpGs overlap with Horvath CpGs
load("/Users/kevinjohnson/Komen_Normal_5mC/05.Genomic_Enrichment/Files/Komen_age_related_CpGs.RData")
sum(Age_DMPs%in%Horvath_CpGs) # 17 out of 353 MEASURED CpGs

############################################################
# Part 5: Compare Epigenetic Clock in Komen vs. TCGA normal adjacent
############################################################
# 97 adjacent-to-tumor normal breast tissue samples were submitted to the Horvath clock
TCGA_Clock_output <- read.csv("06.Epigenetic_Clock/Files/TCGA_NormalAdjacent_Betas_Clock.output.csv", header=TRUE)
# Print summary for Komen and TCGA
summary(covar_clock$AgeAccelerationResidual)
summary(TCGA_Clock_output$AgeAccelerationResidual)
# Examine the difference in AgeAccel between the two tissue types
wilcox.test(covar_clock$AgeAccelerationResidual, TCGA_Clock_output$AgeAccelerationResidual)

