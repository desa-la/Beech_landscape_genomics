######################################################################
# This script is adapted for our data following the code from:
# https://github.com/Capblancq/RDA-genome-scan/blob/master/Script_RDA.R
# Some commented parts are identical to those from Capblancq code
#######################################################################

library(data.table)
library(dplyr)
library(pegas)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(cowplot)
library(corrplot)
library(psych)

imputed <- fread("/home/lazic@ad.vti.bund.de/P1.Genomic.offset/RDA/final.filtered.bim.indep.lfmm_imputed.lfmm")
pop <- read.table("/home/lazic@ad.vti.bund.de/P1.Genomic.offset/RDA/populations.txt", header = T)

#Getting allele frequencies
setDT(imputed)  # Convert the genotypes data frame to a data.table
setDT(pop)           # Convert the population data frame to a data.table
# If the rows align perfectly
imputed[, pop := pop$pop]
# Calculate the mean per SNP column by population, assuming 'pop' is now within imputed.lfmm
allele_freqs <- imputed[, lapply(.SD, function(x) mean(x, na.rm = TRUE) / 2), by = pop]
allele_freqs_nopop <-select(allele_freqs, -pop)

#wza.allele.frequencies (using this file to get the SNP IDs as collumn names)
chr <- read.table("/home/lazic@ad.vti.bund.de/P1.Genomic.offset/RDA/20.frq", header = TRUE, sep = "\t")
chr <- chr %>% mutate(id = paste0(CHROM, "_", POS))
chr_column = chr[, 7, drop = FALSE]
chr_row <- t(chr_column)
chr_row <- as.data.frame(chr_row)
new_column_names <- as.character(unlist(chr_row[1, ]))
colnames(allele_freqs_nopop) <- new_column_names

#Env data
#Input file for environmental data is same as in baypass analysis, extracted from WorldClim database
#setwd("P:/P1_Genomic_offset/1_Intermediate_files/BayPass")
env <- read.table("/home/lazic@ad.vti.bund.de/P1.Genomic.offset/RDA/baypass.environment.data.txt", header = F)
env <- t(env)
Env <- as.data.frame(env)
## Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
#Import population and coordinates
#setwd("P:/P1_Genomic_offset/1_Intermediate_files/RDA")
pop <- read.table("/home/lazic@ad.vti.bund.de/P1.Genomic.offset/RDA/pop.coordinates.txt", header = T)
## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(pop$prov)

#Inferring population structure
## Running a PCA on all genetic markers
pca <- rda(allele_freqs_nopop, scale=T) # PCA in vegan uses the rda() call without any predictors
#Screeplot of the PCA eigenvalues
png("PCA_Eigenvalues_all.snps.png", width = 800, height = 600)
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")
dev.off()
## Neutral population structure table
PCs <- scores(pca, choices=c(1:6), display="sites", scaling=0)
PopStruct <- data.frame(Population = allele_freqs[,1], PCs)
colnames(PopStruct) <- c("Population", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

#selecting variables that are not highly correlated and that are bilogically meaningful
Env <- Env %>% select(1, 2, 3, 5, 6, 7, 8, 9, 12, 15, 18, 19)
#combining
Variables <- data.frame(PopStruct[,-1], Env)
Variables$lat <- pop$lat
Variables$long <- pop$long

##########################################
# Variance partitioning
##########################################
## Full model
pRDAfull <- rda(allele_freqs_nopop ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + long + lat + V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V12 + V15 + V18 + V19, Variables)
RsquareAdj(pRDAfull)
anova(pRDAfull)

## Pure climate model
pRDAclim <- rda(allele_freqs_nopop ~ V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V12 + V15 + V18 + V19 + Condition(long + lat + PC1 + PC2 + PC3 + PC4 + PC5 + PC6), Variables)
RsquareAdj(pRDAclim)
anova(pRDAclim)

## Pure neutral population structure model  
pRDAstruct <- rda(allele_freqs_nopop ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Condition(long + lat + V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V12 + V15 + V18 + V19), Variables)
RsquareAdj(pRDAstruct)
anova(pRDAstruct)

##Pure geography model 
pRDAgeog <- rda(allele_freqs_nopop ~ long + lat + Condition(V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V12 + V15 + V18 + V19 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6), Variables)
RsquareAdj(pRDAgeog)
anova(pRDAgeog)

####################################################################################################
######### Correlation among variables using a correlogram
####################################################################################################
png("corrplot_12variables_PCs_lat_long.png", width = 800, height = 600)
corrplot(cor(Variables[, c("PC1","PC2","PC3","PC4","PC5","PC6","long","lat","V1","V2","V3","V5","V6","V7","V8","V9","V12","V15","V18","V19")]), type="upper")
dev.off()

########################################################################################
##########  Genotype-Environment Associations: identifying loci under selection
########################################################################################

#Conducting the genome scan using pRDA
RDA_env <- rda(allele_freqs_nopop ~ V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V12 + V15 + V18 + V19 + Condition(PC1 + PC2 + PC3 + PC4 + PC5 + PC6),  Variables)
#Conducting the genome scan using simple model without structure correction
RDA_env_nocorrection <- rda(allele_freqs_nopop ~ V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V12 + V15 + V18 + V19,  Variables)

# choose a number of RDA axes to include when conducting the genome scan
png("RDA_Eigenvalues_of_constrained_axes.png", width = 800, height = 600)
screeplot(RDA_env, main="Eigenvalues of constrained axes")
dev.off()
# choose a number of RDA axes to include when conducting the genome scan
png("RDA_Eigenvalues_of_constrained_axes_no_structure_correction.png", width = 800, height = 600)
screeplot(RDA_env_nocorrection, main="Eigenvalues of constrained axes")
dev.off()

#load the rdadapt function, described in Capblancq et al. (2018) and use it to conduct the genome scan
rdadapt<-function(rda,K)
{
	  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
    resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
    lambda <- median(resmaha)/qchisq(0.5,df=K)
      reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
      qval <- qvalue(reschi2test)
        q.values_rdadapt<-qval$qvalues
        return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

## Running the function with K = 2
rdadapt_env<-rdadapt(RDA_env, 2)
rdadapt_env_nocorrection<-rdadapt(RDA_env_nocorrection, 2)
#One critical step when conducting a genome scan is to set a pertinent p-value threshold to identify the outlier loci.
#Here, we used a Bonferroni correction to account for multiple testing.  
#the rdadapt function returns both p-values and q-values, which means it is possible to use a 
#FDR (False Discovery Rate) approach instead of a p-value threshold to identify outliers
## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)
### Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(allele_freqs_nopop)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)])
write.table(outliers, "outliers_rdadapt_env_12variables_structure_corrected.txt", sep = "\t", col.names = NA, quote = FALSE)
#For no structure correction
outliers_no_correction <- data.frame(Loci = colnames(allele_freqs_nopop)[which(rdadapt_env_nocorrection$p.values<thres_env)], p.value = rdadapt_env_nocorrection$p.values[which(rdadapt_env_nocorrection$p.values<thres_env)])
write.table(outliers_no_correction, "outliers_rdadapt_env_12variables_NO_structure_corrected.txt", sep = "\t", col.names = NA, quote = FALSE)

#save.image(file='real_12var.RData')
#load('real_12var.RData')
