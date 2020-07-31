library(tidyverse)
library(bayou)
library(BioGeoBEARS)

par(mfrow = c(1,2))
plotSimmap.mcmc(chain_free, edge.type = "theta", 
                pp.cutoff = .3, cex = .01, no.margin = T, circles = F)
plotSimmap.mcmc(chain_free, edge.type = "regimes", 
                pp.cutoff = .3, cex = .01, no.margin = T, direction = "leftwards")

###
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

plotSimmap.mcmc(chain_free, edge.type = "regimes", 
                pp.cutoff = .3, cex = .01, no.margin = T,  direction = "leftwards")
