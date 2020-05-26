library(GenSA)
library(FD)      
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(tidyverse)

spptraits1 <- read.csv("data/spptraits1.csv")
spptraits1 %>% filter(!is.na(fruit_lg)) %>% select(Gen_spp, sunda, sulawesi, maluku, newguinea) %>% 
  transmute(Gen_spp,loc = paste0(sunda, sulawesi, maluku, newguinea)) %>% 
  write.table(file = "data/indomalay_geog.csv", sep = " ", row.names = F)

# file path
trfn <- "data/zanne_tree.tre"
geog_ind <- "data/indomalay_geog.data"

moref(geog_ind)

tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geog_ind)
max(rowSums(dfnums_to_numeric(tipranges@df))) # Maximum range size observed
max_range_size <- 4 # max numb. areas