setwd("C:/Users/Desanka Lazic/Documents/Beech_Thuenen")
library(dplyr)

beech_pop <- read.table("populations.final.UAcorrected.txt", header = TRUE)
names <- read.table("names_country.txt", header = TRUE)
ordered_samples <- arrange(names, country)
write.table(ordered_samples, file = "populations.txt", sep = " ", quote = FALSE, row.names = FALSE)
provenances <- unique(beech_pop$prov)
provenances <- subset(beech_pop, select=c(genotype, prov))
write.table(provenances, file ="provenances_quotes.txt", quote = TRUE, row.names = FALSE, sep = " ")

#Getting the population file for each provenance (for defining populations in Fst analysis)
provenances <- read.table("names_prov.txt", header = TRUE)
unique_prov <- unique(provenances$prov)

#function to loop through each unique "prov" value and create separate files
for (prov_value in unique_prov) {
  # Subset the data frame for the current "prov" value
  subset_data <- provenances[provenances$prov == prov_value, ]
  
  # Extract genotypes from the subset
  genotypes <- subset_data$genotype_genotype
  
  # Create a filename based on the "prov" value
  filename <- paste0("prov_", prov_value, ".txt")
  
  # Write the genotypes to the file
  write.table(genotypes, file = filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

write.table(unique_prov, file = "uniq_prov.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#After running pairwise Fst across all populations following 2_Fst.sh script in this folder, we proceed with IBD
#####################################################
#################### Mantel #########################
####################################################
library(geosphere)
library(tidyverse)
library(vegan)
library(cowplot)
library(ggplot2)

#import Fst values between pops
setwd("P:/P1_Genomic_offset/IBD_IBE")
fst <- read.table("fst_matrix.txt", header = TRUE)
fst <- fst %>% select(Population1, Population2, Pairwise_fst)

########################################
#create a pairwise matrix
fst_matrix <- fst %>%
  pivot_wider(names_from = Population2, values_from = Pairwise_fst) %>%
  column_to_rownames(var = "Population1")
fst_matrix1 <- data.matrix(fst_matrix)
# Replace NA values with zeros
na_indices <- is.na(fst_matrix1)
fst_matrix1[na_indices] <- 0

#add populations missing due to formating from vcftools
fst_matrix2 <- cbind(pop146 = 0, fst_matrix1)
fst_matrix3 <- rbind(fst_matrix2, pop51 = 0)
# transpose the matrix
fst_matrix4 <- fst_matrix3 + t(fst_matrix3)
#export the matrix
write.table(fst_matrix, file = "fst_matrix_R.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = " ")

#Getting in the coordinates
beech_pop <- read.table("populations.final.UAcorrected.txt", header = TRUE)
#Subset populations/provenances to get coordinate for each pop for the matrix
prov <- distinct(beech_pop, prov, lat, long)
#creating matrix of geographic distances between populations
distance=as.data.frame(prov[,c(3:2)])
rownames(distance)=prov$prov
geo.dist = distm(distance, fun=distVincentyEllipsoid)
rownames(geo.dist) =prov$prov
colnames(geo.dist) =prov$prov
write.csv(as.matrix(geo.dist),file="geo_dist_provenances.csv",quote=F)

##### IBD
mantel(fst_matrix4, geo.dist, method="pearson", permutations=999)

########## Plotting #################
#This function is from: https://github.com/jingwanglab/Populus_genomic_prediction_climate_vulnerability/blob/main/3-Population_genetics/6IBD%26IBE.txt with adapted numbers for our data
######################################
trans <- function(raw_data){
  raw_data=raw_data[-1,-1]
  out_data=data.frame(raw_data[,1])
  colnames(out_data)="value"
  for (i in 1:97){
    temp=data.frame(raw_data[i:97,i])
    colnames(temp)="value"
    out_data=rbind(out_data,temp)
  }
  out_data=na.omit(out_data)
  return(out_data)
}

plot_data=as.data.frame(matrix(nrow=4850,ncol=0))
plot_data$geo.dists=trans(geo.dists)$value
plot_data$fst_matrix4=trans(fst_matrix4)$value
colnames(plot_data)=c("geo_dist","fst")
write.csv(plot_data,file="plot_data.csv",quote=F,row.names = F)

plot_data=read.csv("plot_data.csv",head=T)
plot_data <- plot_data %>% mutate(geo_dist = geo_dist/1000)
plot_data <- plot_data %>% select(geo_dist, fst)

p1=ggplot(plot_data)+
  geom_point(aes(x=geo_dist,y=fst),size = 3,alpha=0.7,color="black",shape=21,fill="#708090")+
  geom_smooth(aes(x=geo_dist, y=fst),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#9c9c9c", fill="#d6d6d6",size = 1.5,fullrange = F) +
  labs(x = "Geographical Distance (km)",y = expression(italic(F)[italic(ST)]),size = 5.5)+
  panel_border(color = "black", size = 0.6, linetype = 1, remove = FALSE)+
  theme_bw()+
  theme(text=element_text(family="sans"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x=element_line(colour="black"),
        #axis.line.y=element_line(colour="black"),
        #panel.border=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) + annotate("text", x=150, y=0.085, label= "Mantels' r = 0.79")
