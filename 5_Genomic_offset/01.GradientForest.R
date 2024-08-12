###################################################################################################################
# This code was partially adapted from github repository:
# https://github.com/jingwanglab/Populus_genomic_prediction_climate_vulnerability/blob/main/7-Genetic_offset/1GF.R
###################################################################################################################

library(gradientForest)
library(LEA)
library(data.table)
library(dplyr)
#Import genetic data
setwd("P:/P1_Genomic_offset/1_Intermediate_files/WZA")
#Imported data via "Import dataset"
#wza.allele.frequencies
chr <- read.table("20.frq", header = TRUE, sep = "\t")
combine <- cbind(chr[, 1:2], wza.allele.frequencies)
# extract frequency data for each SNP
freqs_all <- combine[,3:198]
# 1. choosing only reference allele per population
freqs_all <- freqs_all[, seq(1, ncol(freqs_all), by = 2)]

##########################
tfreqs <- t(freqs_all)
tfreqs <- as.data.frame(tfreqs)
set.seed(123)
# Sample 10,000 row indices
col_indices <- sample(1:ncol(tfreqs), 10000, replace = FALSE)
# Subset the data frame based on the sampled row indices
freq.all_SNPs <- tfreqs[,col_indices]
#Removing populations with less than 5 individuals
rows_to_remove <- c(91, 8, 51, 98, 90, 53)
all.SNPs.noNA <- freq.all_SNPs %>% slice(-rows_to_remove)

#Get the environment data
setwd("P:/P1_Genomic_offset/1_Intermediate_files/GF")
#Choose NON scaled env 
envGF <- read.table("bio.present.env", header = FALSE, sep = " ")

setwd("P:/P3_GWAS/1_Intermediate_files/Phenotyping_related")
UAcorrected <- read.table("populations.final.UAcorrected.txt", sep = " ", header = T)
prov <- UAcorrected[,2]
env <- envGF %>%
  mutate(pop = prov)
pop_based <- env %>% distinct(pop, .keep_all = TRUE)
population_order <- pop_based[,20]
#remove last column
pop_based <- pop_based[, -ncol(pop_based)]
envGF <- pop_based
#Removing the same populations as in genetic data
envGF.noNA <- envGF %>% slice(-rows_to_remove)

#########################################################################
nSpecs <- dim(all.SNPs.noNA)[2]
maxLevel <- log2(0.368*nrow(envGF.noNA)/2)
all_gfmod_maf <- gradientForest(cbind(envGF.noNA, all.SNPs.noNA), predictor.vars=colnames(envGF.noNA),
                                response.vars=colnames(all.SNPs.noNA), ntree=500, maxLevel=maxLevel,
                                trace=T, corr.threshold=0.5)
#########################################################################
setwd("P:/P1_Genomic_offset/3_Results/GF/no_scaled_env")
pdf(file="all_predictoroverallimportance_clean_script.pdf")
plot(all_gfmod_maf,plot.type="O")
dev.off()
#########################################################################

########################################################################################################
#fut45
#Get the environmental data for future climate scenario
setwd("P:/P1_Genomic_offset/1_Intermediate_files/GF")
#NON scaled
env45 <- read.table("bio.fut.45.mpi.env", header = FALSE, sep = " ")
#setwd("P:/P3_GWAS/1_Intermediate_files/Phenotyping_related")
#UAcorrected <- read.table("populations.final.UAcorrected.txt", sep = " ", header = T)
prov <- UAcorrected[,2]
env.45 <- env45 %>%
  mutate(pop = prov)
pop_based_fut <- env.45 %>% distinct(pop, .keep_all = TRUE)
population_order_fut <- pop_based_fut[,20]
#remove last column
pop_based_fut <- pop_based_fut[, -ncol(pop_based_fut)]
env45<- pop_based_fut
#Removing the same populations as in genetic data
env45.noNA <- env45 %>% slice(-rows_to_remove)

#Predict current climate
all_predict=predict(all_gfmod_maf,envGF.noNA[,c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19")])
#Predict future climate
future45_predict=predict(all_gfmod_maf, env45.noNA[,c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19")])
#Calculate Genomic offset using Euclidian distance formula
genOffset45<-sqrt((future45_predict[,1]-all_predict[,1])^2+(future45_predict[,2]-all_predict[,2])^2+(future45_predict[,3]-all_predict[,3])^2+(future45_predict[,4]-all_predict[,4])^2+(future45_predict[,5]-all_predict[,5])^2+(future45_predict[,6]-all_predict[,6])^2+(future45_predict[,7]-all_predict[,7])^2+(future45_predict[,8]-all_predict[,8])^2+(future45_predict[,9]-all_predict[,9])^2+(future45_predict[,10]-all_predict[,10])^2+(future45_predict[,11]-all_predict[,11])^2+(future45_predict[,12]-all_predict[,12])^2+(future45_predict[,13]-all_predict[,13])^2+(future45_predict[,14]-all_predict[,14])^2+(future45_predict[,15]-all_predict[,15])^2+(future45_predict[,16]-all_predict[,16])^2+(future45_predict[,17]-all_predict[,17])^2+(future45_predict[,18]-all_predict[,18])^2+(future45_predict[,19]-all_predict[,19])^2)
setwd("P:/P1_Genomic_offset/3_Results/GF/no_scaled_env")
Offset45=cbind(lat.long,genOffset45)
colnames(Offset45)[4]<-"offset"
write.csv(Offset45, file = "offset_fut45_clean_script.csv", quote=F, row.names=T)
