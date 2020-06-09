library(GenSA)
library(FD)      
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(tidyverse)

# data cleaning  and preparation #####
spptraits1 <- read.csv("data/spptraits1.csv")
geo_tab <- spptraits1 %>% filter(!is.na(fruit_lg)) %>% select(Gen_spp, sunda, sulawesi, maluku, newguinea) %>%
  filter(!is.na(sunda), !is.na(sulawesi), !is.na(maluku), !is.na(newguinea)) %>% 
  transmute(Gen_spp, loc = paste0(sunda, sulawesi, maluku, newguinea)) %>% filter(loc != "0000")

spptree1 <- read.tree("data/zanne_tree.tre")
pr_tree <- drop.tip(spptree1, setdiff(spptree1$tip.label, geo_tab$Gen_spp))
pr_tree$tip.label <- gsub("-", "_", pr_tree$tip.label)
write.tree(pr_tree,"data/zanne_tree_pr.tre")

geo_tab$Gen_spp <- gsub("-", "_", geo_tab$Gen_spp)
geo_data <- geo_tab[!is.na(match(geo_tab$Gen_spp, pr_tree$tip.label)),]

setdiff(geo_data$Gen_spp, pr_tree$tip.label)
setdiff(pr_tree$tip.label, geo_data$Gen_spp)
duplicated(geo_data$Gen_spp) %>% table()
duplicated(pr_tree$tip.label) %>% table()

write.table(geo_data, file = "data/indomalay_geog.data", sep = "\t", row.names = F)

# file path
trfn <- "data/zanne_tree_pr.tre"
geogfn <- "data/indomalay_geog.data"

moref(trfn)
moref(geogfn)

tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geog_ind)
max(rowSums(dfnums_to_numeric(tipranges@df))) # Maximum range size observed
max_range_size <- 4 # max numb. areas

# number of states
numstates_from_numareas(numareas = 4, maxareas = 4, include_null_range = TRUE)

# DEC ####
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size

BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc

# files to use in time-stratified analyses
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"

BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE   # force_sparse=TRUE causes pathology & isn't much faster at this scale

# dispersal multiplier matrix load
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
#BioGeoBEARS_run_object$master_table # stratified tree described in table

# Default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run

check_BioGeoBEARS_run(BioGeoBEARS_run_object) # check for errors

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow <- TRUE
resfn <- "imdomalay_DEC_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res <- bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file = resfn)
  resDEC <- res
} else {
  # Loads to "res"
  load(resfn)
  resDEC <- res
}

# DEC + J ####
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    
BioGeoBEARS_run_object$include_null_range <- TRUE   

BioGeoBEARS_run_object$on_NaN_error <- -1e50    
BioGeoBEARS_run_object$speedup <- TRUE         
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE 

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE   

# Set up DEC+J model, get ML par vals  2-parameter nested model
dstart <- resDEC$outputs@params_table["d","est"]
estart <- resDEC$outputs@params_table["e","est"]
jstart <- 0.0001

# d, e starting vals
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"]<- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"]<- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] <- estart

# j as free par
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]  <- jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn <- "imdomalay_DEC+J_M0_unconstrained_v1.Rdata"
runslow <- TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  res <- bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file = resfn)
  
  resDECj <- res
} else {
  # Loads to "res"
  load(resfn)
  resDECj <- res
}

# DEC and DEC+J plots ####
pdffn <- "indomalay_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width = 6, height = 6)

analysis_titletxt <- "indomalay DEC M0_unconstrained"
results_object <- resDEC
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
res2 <- plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                 addl_params = list("j"), plotwhat = "text", 
                                 label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                                 plotsplits = TRUE, cornercoords_loc = scriptdir, 
                                 include_null_range = TRUE, tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                         plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tr, tipranges = tipranges)


analysis_titletxt <- "indomalay DEC+J M0_unconstrained"
results_object <- resDECj
res1 <- plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                 addl_params = list("j"), plotwhat = "text", 
                                 label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                                 plotsplits = TRUE, cornercoords_loc = scriptdir, 
                                 include_null_range = TRUE, tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                         plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tr, tipranges = tipranges)
dev.off()
system(paste("open ", pdffn, sep = ""))

# DIVA-like ####
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE  

BioGeoBEARS_run_object$on_NaN_error <- -1e50    
BioGeoBEARS_run_object$speedup <- TRUE         
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE 

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE   

# setup DIVA-like model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "ysv*1/2"

# widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"]  <- 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow <- TRUE
resfn <- "indomalay_DIVALIKE_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res <- bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file = resfn)
  resDIVALIKE <- res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE <- res
}

# DIVA-like+J ####
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE  

BioGeoBEARS_run_object$on_NaN_error <- -1e50    
BioGeoBEARS_run_object$speedup <- TRUE         
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE 

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE   

# Set up DIVALIKE+J model
dstart <- resDIVALIKE$outputs@params_table["d","est"]
estart <- resDIVALIKE$outputs@params_table["e","est"]
jstart <- 0.0001

# d, e starting values 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"]  <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"]  <- estart

# no subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"]  <- 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"]  <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"]   <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"]   <- "ysv*1/2"

# widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"]  <- 0.5

# jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]  <- jstart

# "j" max should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] <- 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)
resfn <- "indomalay_DIVALIKE+J_M0_unconstrained_v1.Rdata"
runslow <- TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res <- bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file = resfn)
  
  resDIVALIKEj <- res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj <- res
}

# DIVA-like and DIVA-like+J plots ####
pdffn <- "indomalay_DIVALIKE_vs_DIVALIKE+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width = 6, height = 6)

analysis_titletxt <- "DIVALIKE indomalay M0_unconstrained"
results_object <- resDIVALIKE
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))

res2 <- plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                 addl_params = list("j"), plotwhat = "text", 
                                 label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                                 plotsplits = TRUE, cornercoords_loc = scriptdir, 
                                 include_null_range = TRUE, tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                         plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tr, tipranges = tipranges)

analysis_titletxt <- "DIVALIKE+J indomalay M0_unconstrained"
results_object <- resDIVALIKEj

res1 <- plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                 addl_params = list("j"), plotwhat = "text", 
                                 label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                                 plotsplits = TRUE, cornercoords_loc = scriptdir, 
                                 include_null_range = TRUE, tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.7, statecex = 0.7, splitcex = 0.6, titlecex = 0.8, 
                         plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tr, tipranges = tipranges)
dev.off()
system(paste("open ", pdffn, sep = ""))

# BAYEAREA-like ####
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE  

BioGeoBEARS_run_object$on_NaN_error <- -1e50    
BioGeoBEARS_run_object$speedup <- TRUE         
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE 

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE 

# setup BAYEAREA-like model
# no subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"]  <- 0.0

# no vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"]  <- 0.0

# no jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]  <- 0.01

# linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"]  <- "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"]   <- "1-j"

# Only sympatric/range-copying (y) events allowed, exact copying (both descendants , same size as ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"]  <- 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)
runslow <- TRUE
resfn <- "indomalay_BAYAREALIKE_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res <- bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file = resfn)
  resBAYAREALIKE <- res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE <- res
}

# BAYEAREA-like+J ####
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE  

BioGeoBEARS_run_object$on_NaN_error <- -1e50    
BioGeoBEARS_run_object$speedup <- TRUE         
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE 

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE 

# set up BAYAREA-like+J model
dstart <- resBAYAREALIKE$outputs@params_table["d","est"]
estart <- resBAYAREALIKE$outputs@params_table["e","est"]
jstart <- 0.0001

# d, e starting values 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"]  <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"]  <- estart

# no subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"]  <- 0.0

# no vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"]  <- 0.0

#  jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]  <- jstart

# BAYAREA-like+J, "j" max should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 0.99999

# linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"]  <- "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"]   <- "1-j"

# Only sympatric/range-copying (y) events allowed, exact copying (both descendants , same size as ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"]  <- 0.9999

##### if computer crashes !!###
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn <- "indomalay_BAYAREALIKE+J_M0_unconstrained_v1.Rdata"
runslow <- TRUE
if (runslow)
{
  res <- bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file = resfn)
  
  resBAYAREALIKEj <- res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj <- res
}

# sum stats ####

restable <- NULL
teststable <- NULL

## DEC & DEC+J ##
LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 <- 3 # this can be fixed putting independet names to results given the model !!!
numparams2 <- 2
stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC is the null model
res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC, 
                                                       returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECj, 
                                                       returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)

tmp_tests <- conditional_format_table(stats)

restable <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)

## DIVA-like & DIVA-like+J ##
LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 <- 3
numparams2 <- 2
stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKE, 
                                                       returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKEj, 
                                                       returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
tmp_tests <- conditional_format_table(stats)
restable <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)

## BAYEAREA-like & BAYEAREA-like+J ##
LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 <- 3
numparams2 <- 2
stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# BAYEAREA-like, null model for Likelihood Ratio Test (LRT)
res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKE, 
                                                       returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
# BAYEAREA-like+J, alternative model for Likelihood Ratio Test (LRT)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKEj, 
                                                       returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
tmp_tests <- conditional_format_table(stats)
restable <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)

# asseembly
teststable$alt <- c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null <- c("DEC", "DIVALIKE", "BAYAREALIKE")

#do.call(cbind, teststable) %>% str()
#t_temp <- map_df(teststable, ~as.data.frame(t(.))) %>% str() %>% data.frame()
#colnames(t_temp) <- colnames(teststable)
#unlist(t_temp)

row.names(restable) <- c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable <- put_jcol_after_ecol(restable)

write.table(restable, "results_table.csv", sep = "\t", row.names = F)
#write.table(teststable, "test_table.csv", sep = "\t", row.names = F)


## model weigths ##
restable2 <- restable

# AICs
AICtable <- calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable <- cbind(restable, AICtable)
restable_AIC_rellike <- AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike <- put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike %>% round(3) 
write.table(restable_AIC_rellike, "results_table_relLik.csv", sep = "\t", row.names = F)

# AICcs
samplesize <- length(tr$tip.label)
AICtable <- calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 <- cbind(restable2, AICtable)
restable_AICc_rellike <- AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike <- put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike %>% round(3)
write.table(restable_AICc_rellike, "results_table_rel_Lik_cor.csv", sep = "\t", row.names = F)
