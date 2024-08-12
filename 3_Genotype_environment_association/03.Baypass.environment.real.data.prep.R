library(tidyverse)
library(dplyr)
library(ggplot2)

#Setting up environment data for baypass
#Upload the file
setwd("P:/P3_GWAS/1_Intermediate_files/Phenotyping_related")
#Get the envirnonmental data extracted from WorldClim database
env <- read.table("bio.present.env", sep = " ")
#Get population info
UAcorrected <- read.table("populations.final.UAcorrected.txt", sep = " ", header = T)
prov <- UAcorrected[,2]
env <- env %>%
  mutate(pop = prov)
pop_based <- env %>% distinct(pop, .keep_all = TRUE)
population_order <- pop_based[,20]
#remove last column
pop_based <- pop_based[, -ncol(pop_based)]
#Transposing data frame
library(data.table)
pop.based.t <- transpose(as.data.table(pop_based))
pop.based.t <- t(pop_based)
setwd("P:/P1_Genomic_offset/1_Intermediate_files/BayPass")
write.table(pop.based.t,
            file = "baypass.environment.data.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE, quote = FALSE)
