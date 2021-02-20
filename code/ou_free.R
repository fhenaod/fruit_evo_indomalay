library(dplyr)
library(treeplyr)
library(bayou)

# data
spptree1<-read.tree("data/zanne_tree.tre")
spptraits1<-read.csv("data/spptraits1.csv")

tree_data <- make.treedata(spptree1, spptraits1)
tree_data_ln <- mutate(tree_data, ln_fruit_lg = log(fruit_lg), ln_seed_lg = log(seed_lg), ln_num_seeds = log(num_seeds))
sd2 <- sqrt(log(1+ (var(pull(tree_data_ln$dat,fruit_lg), na.rm = T)) / (mean(pull(tree_data_ln$dat,fruit_lg), na.rm = T))^2 ))

d_fruit_lg_ln <- (filter(tree_data_ln, !is.na(ln_fruit_lg)))
name.check(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg))

# analys
priorOU <- make.prior(d_fruit_lg_ln$phy, 
                    dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dk = "cdpois", dtheta = "dnorm"),
                    param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                 dsb = list(bmax = 1, prob = 1), 
                                 dk = list(lambda = (2*Ntip(d_fruit_lg_ln$phy)-2)*(2.5/100), kmax = (2*Ntip(d_fruit_lg_ln$phy)-2)*(5/100)),
                                 dtheta = list(mean = mean(getVector(d_fruit_lg_ln, ln_fruit_lg)), 
                                               sd = sd2))
)

mcmcOU<-bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg),
                       prior = priorOU, new.dir = "modelOU/", outname = "modelOU_r001", plot.freq = NULL) 

gens<-2000
mcmcOU$run(gens)
saveRDS(mcmcOU, file = "modelOU/mcmcOU.rds")
shiftsum <- shiftSummaries(chain, mcmcOU)
plotShiftSummaries(shiftsum)

setwd("chain2/")
mcmc <- readRDS("model_free/mcmcOU.rds")

# PP >= 0.3
chain_free <- mcmc$load()
chain_free <- set.burnin(chain_free, .3)
saveRDS(chain_free, file = "model_free/chain_free_postBI.rds")
chain_free <- readRDS("model_free/chain_free_postBI.rds")
plotSimmap.mcmc(chain_free, pp.cutoff = .3, cex = .01, no.margin = T)

sum_cf <- summary(chain_free)
saveRDS(sum_cf, file = "model_free/chain_sum.rds")
sum_cf <- readRDS("model_free/chain_sum.rds")

shiftsum <- shiftSummaries(chain_free, mcmc, pp.cutoff = .3)
saveRDS(shiftsum, file = "model_free/shiftsum_free.rds")
shiftsum <- readRDS("chain2/model_free/shiftsum_free.rds")
pdf(paste0("shiftsummaryplot.pdf"))
plotShiftSummaries((shiftsum))
dev.off()

# PP >= 0.5 & 0.75
mcmcOU <- readRDS("ou_free/chain2_pp.5/model_free/mcmcOU.rds")
chainOU <- readRDS("ou_free/chain2_pp.5/model_free/mcmcOU_chain.rds")

sum_cf <- readRDS("ou_free/chain2_pp.5/model_free/chain_sum.rds")

shiftsum.5 <- readRDS("ou_free/chain2_pp.5/model_free/shiftsum_free.5.rds")
shiftsum.75 <- readRDS("ou_free/chain2_pp.5/model_free/shiftsum_free.75.rds")

# to check OU model adequacy ####
ou_fit<- fitContinuous(shiftsum$tree, shiftsum$dat, model = "OU", SE = sd/10)
#cont<-ape::pic(shiftsum$dat, shiftsum$tree, scaled = T)
arb_fit<-arbutus::arbutus(ou_fit)
arb_fit
plot(arb_fit)
