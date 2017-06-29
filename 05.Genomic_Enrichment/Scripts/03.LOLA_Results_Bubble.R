#####################################################
# DNA methylation differences at regulatory elements are associated with the cancer risk factor age in normal breast tissue 
# Author: Kevin Johnson, Andy Houseman, Jessica King, Brock Christensen
######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Bubble graph for LOLA TFBS enrichment results (ENCODE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
# Use ggrepel package to distance text away from other points
# devtools::install_github("slowkow/ggrepel")
library(ggrepel)

# Project working directory
setwd("/Users/kevinjohnson/Komen_Normal_5mC/")

#######################
# Age-related hypermethylated
#######################
LOLA_TFBS_hyper = read.table("05.Genomic_Enrichment/lolaResults_hyper/col_encode_tfbs.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
# Select mammary relevant cell-types
LOLA_TFBS_hyper_sub <- LOLA_TFBS_hyper[LOLA_TFBS_hyper$cellType=='HMEC' | LOLA_TFBS_hyper$cellType=='MCF-7' | LOLA_TFBS_hyper$cellType=='MCF10A-Er-Src' | LOLA_TFBS_hyper$cellType=='T-47D',]
LOLA_TFBS_hyper_sub$TF <- unlist(lapply(strsplit(LOLA_TFBS_hyper_sub$antibody, "_"), '[[', 1))
LOLA_TFBS_hyper_sub$TF <- as.factor(LOLA_TFBS_hyper_sub$TF)
nom <- -log10(0.01)
-log10(LOLA_TFBS_hyper_sub$qValue)
# Off-set the q-values of 1.
LOLA_TFBS_hyper_sub$qValue_adj = ifelse(LOLA_TFBS_hyper_sub$qValue==1.0, LOLA_TFBS_hyper_sub$qValue-0.00001, LOLA_TFBS_hyper_sub$qValue)

png("05.Genomic_Enrichment/Figures/TFBS_hyper_LOLA.png", height=7*300, width=7*300, res=300)
overall_theme <-  theme(axis.text=element_text(size=12, face="bold"),
                        axis.title=element_text(size=14, face="bold"), axis.text.x = element_text(angle = 45, hjust = 1),
                        panel.grid.major = element_line(colour = "grey40"),
                        panel.grid.minor = element_blank(), legend.key = element_blank())
 
TFBS <- ggplot(LOLA_TFBS_hyper_sub, aes(x = TF, y = -log10(qValue), label = TF)) +
  geom_point(aes(size = logOddsRatio, colour = cellType, alpha=.5)) + 
  geom_text_repel(size = 3) + ylim(0,30) +
  scale_size(range = c(1,15))  +
  geom_hline(yintercept = nom, color = "red", linetype = "dashed", size = 1.1) 
TFBS + theme_bw() + overall_theme + scale_alpha(guide = 'none') + xlab("Transcription Factor") + ylab("-log10(Q-value)")
dev.off()

# How many ChIP-seq experiments with a Q-value < 0.01?
sum(LOLA_TFBS_hyper_sub$qValue<0.01) # 14 different experimental conditions


#######################
# Age-related hypomethylated
#######################
LOLA_TFBS_hypo = read.table("05.Genomic_Enrichment/lolaResults_hypo/col_encode_tfbs.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
# Select mammary associated ENCODE cell lines
LOLA_TFBS_hypo_sub <- LOLA_TFBS_hypo[LOLA_TFBS_hypo$cellType=='HMEC' | LOLA_TFBS_hypo$cellType=='MCF-7' | LOLA_TFBS_hypo$cellType=='MCF10A-Er-Src' | LOLA_TFBS_hypo$cellType=='T-47D',]
LOLA_TFBS_hypo_sub$TF <- unlist(lapply(strsplit(LOLA_TFBS_hypo_sub$antibody, "_"), '[[', 1))
LOLA_TFBS_hypo_sub$TF <- as.factor(LOLA_TFBS_hypo_sub$TF)
LOLA_TFBS_hypo_sub$cellType <- as.factor(LOLA_TFBS_hypo_sub$cellType)
nom <- -log10(0.01)

png("05.Genomic_Enrichment/Figures/TFBS_hypo_LOLA.png", height=7*300, width=7*300, res=300)
overall_theme <-  theme(axis.text=element_text(size=12, face="bold"),
                        axis.title=element_text(size=14, face="bold"), axis.text.x = element_text(angle = 45, hjust = 1),
                        panel.grid.major = element_line(colour = "grey40"),
                        panel.grid.minor = element_blank(), legend.key = element_blank())

TFBS <- ggplot(LOLA_TFBS_hypo_sub, aes(x = TF, y = -log10(qValue),label = TF)) +
  geom_point(aes(size = logOddsRatio, colour = cellType, alpha=.5)) + 
  geom_text_repel(size = 3) + ylim(0,30) +
  scale_size(range = c(1,15))  +
  geom_hline(yintercept = nom, color = "red", linetype = "dashed", size = 1.1) 
TFBS + theme_bw() + overall_theme + scale_alpha(guide = 'none') + xlab("Transcription Factor") + ylab("-log10(Q-value)")
dev.off()

# How many ChIP-seq experiments with a Q-value < 0.01?
sum(LOLA_TFBS_hypo_sub$qValue<0.01) # 8 different experimental conditions


# Write out subseted hyper and hypo enrichment sets for Supplemental Table
write.csv(LOLA_TFBS_hyper_sub, "05.Genomic_Enrichment/Files/hyper-tfbs-lola.csv")
write.csv(LOLA_TFBS_hypo_sub, "05.Genomic_Enrichment/Files/hypo-tfbs-lola.csv")

