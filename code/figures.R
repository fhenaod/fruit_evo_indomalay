library(tidyverse)
library(bayou)
library(BioGeoBEARS)
library(castor)

chain_fixed <- readRDS("ou_free_fix/model_free/mcmc_fixed.pp75.rds")

# Fig 2 ####
data_spp <- read.csv(file = "data/species_clades.csv", head = T)

par(mfrow = c(1,2))
plotSimmap.mcmc(chain_fixed, edge.type = "theta",
                pp.cutoff = .3, cex = .01, no.margin = T, circles = F,
                pal = terrain.colors, 
                legend_settings = 
                  list(plot = T, x = 0.04, y = 1200, 
                       cex.lab = .5, adjx = .02))

tree2paint <- d_fruit_lg_ln$phy
for(i in 1:length(sup_clad_nod)){
  clad2name <- data_spp %>% filter(superclade == sup_clad_nod[i]) %>% 
    select(species) %>% c()
  
  tree2paint <- paintSubTree(tree2paint, 
               get_mrca_of_set(tree2paint, clad2name$species), 
               state = sup_clad_nod[i], anc.state = "0", stem = T)
}

cols <- setNames(c("grey", viridis::viridis(length(sup_clad_nod))), c("0", sup_clad_nod))
plot(tree2paint, fsize = 0.001, 
     colors = cols, direction = "leftwards")

add.simmap.legend(colors = setNames(cols[2:length(cols)], sup_clad_nod), 
                  fsize = 0.35, prompt = FALSE, x = .8, y = 1540)  

# if tree in fan type
par(mfrow = c(1, 1))
plotSimmap.mcmc(chain_free, edge.type = "theta", type = "fan",
                pp.cutoff = .3, cex = .01, no.margin = T, circles = F,
                pal = terrain.colors, 
                legend_settings = 
                  list(plot = T, x = -1.1, y = .4, 
                       cex.lab = .5, adjx = .02))

sup_clad_nod <- unique(data_spp$superclade)
# clade names in nodes 
for(i in 1:length(sup_clad_nod)){
  clad2name <- data_spp %>% filter(superclade == sup_clad_nod[i]) %>% 
    select(species) %>% c()
  
  nodelabels(sup_clad_nod[i], get_mrca_of_set(d_fruit_lg_ln$phy, clad2name$species),
             #adj = c(1.1, 0),
             frame = "n", cex = 1, font = 3)
}  

# clade names in lines arround tips
cols <- setNames(c(viridis::viridis(length(sup_clad_nod))), c(sup_clad_nod))
cols <- setNames(c(RColorBrewer::brewer.pal(n = length(sup_clad_nod), name = "Set1")), c(sup_clad_nod))

for(i in 1:length(sup_clad_nod)){
  clad2name <- data_spp %>% filter(superclade == sup_clad_nod[i]) %>% 
    select(species) %>% c()
  
  arc.cladelabels(d_fruit_lg_ln$phy, sup_clad_nod[i], get_mrca_of_set(d_fruit_lg_ln$phy, clad2name$species),
             orientation = "horizontal", mark.node = T,
             ln.offset = 1.02, lab.offset = 1.04, cex = .9, lwd = 4,
             wing.length = 0.001, col = cols[i])
} 

# Fig 4 ####
par(mfrow = c(1,2))
analysis_titletxt <- ""
results_object <- resDECj

plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie",
                         tipcex = 0, 
                         label.offset = 0.45, show.tip.label = F, 
                         statecex = 0.4, splitcex = 0.6, titlecex = 0.001,
                         plotsplits = F, cornercoords_loc = scriptdir, 
                         include_null_range = T, tr = tr, tipranges = tipranges, 
                         plotlegend = F, tipboxes_TF = F)

plotSimmap.mcmc(chain_fixed, edge.type = "regimes", 
                pp.cutoff = 0.9, cex = .01, no.margin = T, 
                direction = "leftwards", circles = T)


