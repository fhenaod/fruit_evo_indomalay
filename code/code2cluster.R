library(dplyr)
library(treeplyr)
library(bayou)

spptree1<-read.tree("data/zanne_tree.tre")
spptraits1<-read.csv("data/spptraits1.csv")

tree_data<-make.treedata(spptree1, spptraits1)
tree_data_ln<-mutate(tree_data, ln_fruit_lg = log(fruit_lg), ln_seed_lg = log(seed_lg), ln_num_seeds = log(num_seeds))
sd2<-sqrt(log(1+ (var(pull(tree_data_ln$dat,fruit_lg), na.rm = T)) / (mean(pull(tree_data_ln$dat,fruit_lg), na.rm = T))^2 ))

d_fruit_lg_ln<-(filter(tree_data_ln, !is.na(ln_fruit_lg)))
name.check(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg))

priorOU<-make.prior(d_fruit_lg_ln$phy, 
                    dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dk = "cdpois", dtheta = "dnorm"),
                    param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                 dsb = list(bmax = 1, prob = 1), 
                                 dk = list(lambda = (2*Ntip(d_fruit_lg_ln$phy)-2)*(2.5/100), kmax = (2*Ntip(d_fruit_lg_ln$phy)-2)*(5/100)),
                                 dtheta = list(mean = mean(getVector(d_fruit_lg_ln, ln_fruit_lg)), 
                                               sd = sd2))
)

mcmcOU<-bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg),
                       prior = priorOU, new.dir = "modelOU/", outname = "modelOU_r001", plot.freq = NULL) 

saveRDS(mcmcOU, file = "modelOU/mcmcOU.rds")
mcmcOU$run(2000)
chainOU<-mcmcOU$load(saveRDS = T, file = 'modelOU/mcmcOU_chain.rds')
chainOU<-set.burnin(chainOU, 0.3)
summary(chainOU)
