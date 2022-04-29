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
tr <- read.tree("data/zanne_tree_pr.tre")
geogfn <- "data/indomalay_geog.data"

moref(trfn)
moref(geogfn)

tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
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
runslow <- F
resfn <- "biogeob/imdomalay_DEC_M0_unconstrained_v1.Rdata"
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

resfn <- "biogeob/imdomalay_DEC+J_M0_unconstrained_v1.Rdata"
runslow <- F
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
tr <- read.tree("data/zanne_tree_pr.tre")

pdffn <- "indomalay_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width = 6, height = 6)

analysis_titletxt <- "indomalay DEC M0_unconstrained"
results_object <- resDEC
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
res2 <- plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                 addl_params = list("j"), plotwhat = "text", 
                                 label.offset = 0.45, tipcex = 0.2, statecex = 0.4, 
                                 splitcex = 0.4, titlecex = 0.8, 
                                 plotsplits = TRUE, cornercoords_loc = scriptdir, 
                                 include_null_range = TRUE, tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.2, statecex = 0.4, 
                         splitcex = 0.4, titlecex = 0.8, 
                         plotsplits = F, plotlegend = T, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE)

analysis_titletxt <- "indomalay DEC+J M0_unconstrained"
results_object <- resDECj
res1 <- plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                 addl_params = list("j"), plotwhat = "text", 
                                 label.offset = 0.45, tipcex = 0.4, 
                                 statecex = 0.4, splitcex = 0.6, titlecex = 0.8, 
                                 plotsplits = F, cornercoords_loc = scriptdir, 
                                 include_null_range = TRUE, tr = tr, tipranges = tipranges)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie", 
                         label.offset = 0.45, show.tip.label = F, 
                         statecex = 0.4, splitcex = 0.6, titlecex = 0.7, 
                         plotsplits = F, cornercoords_loc = scriptdir, 
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

runslow <- F
resfn <- "biogeob/indomalay_DIVALIKE_M0_unconstrained_v1.Rdata"
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
resfn <- "biogeob/indomalay_DIVALIKE+J_M0_unconstrained_v1.Rdata"
runslow <- F
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
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- 0.01
#BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"]  <- 0.01

# linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"]  <- "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"]   <- "1-j"

# Only sympatric/range-copying (y) events allowed, exact copying (both descendants , same size as ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"]  <- 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)
runslow <- F
resfn <- "biogeob/indomalay_BAYAREALIKE_M0_unconstrained_v1.Rdata"
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

resfn <- "biogeob/indomalay_BAYAREALIKE+J_M0_unconstrained_v1.Rdata"
runslow <- F
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

# BAYEAREA-like and BAYEAREA-like+J plots #####
pdffn <- "indomalay_BAYAREALIKE_vs_BAYAREALIKEE+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width = 6, height = 6)

analysis_titletxt <- "BAYAREALIKE indomalay M0_unconstrained"
results_object <- resBAYAREALIKE
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

analysis_titletxt <- "BAYAREALIKE+J indomalay M0_unconstrained"
results_object <- resBAYAREALIKEj

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

# Ancestral Range estimation sum stats ####
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

write.table(restable, "biogeob/results_table.csv", sep = "\t")
#write.table(teststable, "test_table.csv", sep = "\t")


## model weigths ##
restable2 <- restable

# AICs
AICtable <- calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable <- cbind(restable, AICtable)
restable_AIC_rellike <- AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike <- put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike %>% round(4) %>% mutate(model = rownames(restable_AIC_rellike)) %>% 
  arrange(desc(AIC_wt)) %>% write.table("biogeob/results_table_relLik.csv", sep = "\t")

# AICcs
samplesize <- length(tr$tip.label)
AICtable <- calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 <- cbind(restable2, AICtable)
restable_AICc_rellike <- AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike <- put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike %>% round(4) %>% mutate(model = rownames(restable_AICc_rellike)) %>% 
  arrange(desc(AICc_wt)) %>% write.table("biogeob/results_table_rel_Lik_cor.csv", sep = "\t")


# non-stratified BSM ####
model_name = "resDECj"
res = resDECj

results_object = res
scriptdir = np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))

res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                                addl_params = list("j"), plotwhat = "text", 
                                label.offset = 0.45, tipcex = 0.7, statecex = 0.3, 
                                splitcex = 0.4, titlecex = 0.5, plotsplits = TRUE, 
                                cornercoords_loc = scriptdir, include_null_range = TRUE, 
                                tr = tr, tipranges = tipranges)

clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

# input data
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = T

if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res = res)
  save(stochastic_mapping_inputs_list, file = BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs 
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed = as.numeric(Sys.time()))

# run analysis
runBSMslow = F
if (runBSMslow == TRUE)
{
  BSM_output = runBSM(res, 
                      stochastic_mapping_inputs_list = stochastic_mapping_inputs_list, 
                      maxnum_maps_to_try = 100, nummaps_goal = 10, 
                      maxtries_per_branch = 100, save_after_every_try = TRUE, 
                      savedir = getwd(), seedval = 12345, 
                      wait_before_save = 0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  load(file = "biogeob/bsm/RES_clado_events_tables.Rdata")
  load(file = "biogeob/bsm/RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# BSM output
clado_events_tables <- BSM_output$RES_clado_events_tables
ana_events_tables <- BSM_output$RES_ana_events_tables

head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames

states_list_0based <- rcpp_areas_list_to_states_list(areas = areas, 
                                                    maxareas = max_range_size, 
                                                    include_null_range = include_null_range)

colors_list_for_states <- get_colors_for_states_list_0based(areanames = areanames, 
                                                           states_list_0based = states_list_0based, 
                                                           max_range_size = max_range_size, plot_null_range = TRUE)
# plot a SM
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
stratified <- FALSE
clado_events_table <- clado_events_tables[[1]]
ana_events_table <- ana_events_tables[[1]]

cols_to_get <- names(clado_events_table[,-ncol(clado_events_table)])
colnums <- match(cols_to_get, names(ana_events_table))
ana_events_table_cols_to_add <- ana_events_table[,colnums]
anagenetic_events_txt_below_node <- rep("none", nrow(ana_events_table_cols_to_add))
ana_events_table_cols_to_add <- cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
rows_to_get_TF <- ana_events_table_cols_to_add$node <= length(tr$tip.label)
master_table_cladogenetic_events <- rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

# if in PDF
#pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
#pdf(file = pdffn, width = 6, height = 6)

# BSM to res objt
master_table_cladogenetic_events <- clado_events_tables[[1]]
resmod <- stochastic_map_states_into_res(res = res, 
                                        master_table_cladogenetic_events = master_table_cladogenetic_events, 
                                        stratified = stratified)

plot_BioGeoBEARS_results(results_object = resmod, 
                         analysis_titletxt = "BSM", 
                         addl_params = list("j"), #label.offset = 0.5, 
                         plotwhat = "pie", cornercoords_loc = scriptdir, 
                         root.edge = TRUE, colors_list_for_states = colors_list_for_states, 
                         skiptree = FALSE, show.tip.label = F)

# Paint on the branch states
paint_stochastic_map_branches(res = resmod, 
                              master_table_cladogenetic_events = master_table_cladogenetic_events, 
                              colors_list_for_states = colors_list_for_states, 
                              lwd = 5, lty = par("lty"), root.edge = TRUE, 
                              stratified = stratified)

plot_BioGeoBEARS_results(results_object = resmod, 
                         analysis_titletxt = "Stochastic map", 
                         addl_params = list("j"), plotwhat = "text", 
                         cornercoords_loc = scriptdir, root.edge = TRUE, 
                         colors_list_for_states = colors_list_for_states, 
                         skiptree = TRUE, show.tip.label = F)
 
# close pdf
#dev.off()
#cmdstr = paste("open ", pdffn, sep = "")
#system(cmdstr)

# summarize SM tables
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames <- names(tipranges@df)
actual_names <- areanames

# Get the dmat and times (if any)
dmat_times <- get_dmat_times_from_res(res = res, numstates = NULL)
dmat_times <- get_dmat_times_from_res(res = res, 
                                     numstates = numstates_from_numareas(numareas = 4, maxareas = 4, include_null_range = TRUE))
dmat_times

# Extract BSM output
clado_events_tables <- BSM_output$RES_clado_events_tables
ana_events_tables <- BSM_output$RES_ana_events_tables

# sim the source areas
BSMs_w_sourceAreas <- simulate_source_areas_ana_clado(res, 
                                                     clado_events_tables, 
                                                     ana_events_tables, 
                                                     areanames)

clado_events_tables <- BSMs_w_sourceAreas$clado_events_tables
ana_events_tables <- BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list <- count_ana_clado_events(clado_events_tables, 
                                     ana_events_tables, areanames, actual_names)

summary_counts_BSMs <- counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Hist event counts
hist_event_counts(counts_list, 
                  pdffn = paste0(model_name, "_histograms_of_event_counts.pdf"))

# to check if ML ancestral state/range probs and mean BSMs add up
library(MultinomialCI) # 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, model_name, 
                tr = NULL, plot_each_node = FALSE, linreg_plot = TRUE, 
                MultinomialCI = TRUE)
