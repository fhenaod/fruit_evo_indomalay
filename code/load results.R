library(dplyr)
library(treeplyr)
library(bayou)

setwd("chain1/")
mcmcOU<-readRDS("modelOU/mcmcOU.rds")
chainOU<-mcmcOU$load()

#or

chainOU<-readRDS("chain1/modelOU/mcmcOU_chain.rds")
chainOU<-set.burnin(chainOU, 0.3)
summary(chainOU)
plot(chainOU, auto.layout = FALSE)

plotSimmap.mcmc(chainOU, burnin = 0.3, pp.cutoff = 0.3, cex = .1, type = "fan")
plotBranchHeatMap(d_fruit_lg_ln$phy, chainOU, "theta", burnin = 0.3, pal = cm.colors, cex = .1, type = "fan")
phenogram.density(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), burnin = 0.3, chainOU, pp.cutoff = 0.3)
