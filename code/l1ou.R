library(dplyr)
library(treeplyr)
library(genlasso)
library(l1ou)

spptree1<-read.tree("data/zanne_tree.tre")
spptraits1<-read.csv("data/spptraits1.csv")

tree_data<-make.treedata(spptree1, spptraits1)
tree_data_ln<-mutate(tree_data, ln_fruit_lg = log(fruit_lg), ln_seed_lg = log(seed_lg), ln_num_seeds = log(num_seeds))

d_fruit_lg_ln<-(filter(tree_data_ln, !is.na(ln_fruit_lg)))
name.check(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg))

fruit_lg_adj<-adjust_data(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg))

eModel<-estimate_shift_configuration(fruit_lg_adj$tree, fruit_lg_adj$Y, nCores = 10)
configuration_ic(lizard$tree, eModel$Y, eModel$shift.configuration, criterion = "pBIC") # phylogenetic BIC

stopifnot( identical(sort(eModel$shift.configuration), sort(eModel$shift.configuration)))
saveRDS(eModel, "output/eModel.rds")

nEdges<-Nedge(fruit_lg_adj$tree)
ew<-rep(1, nEdges)
ew[eModel$shift.configuration]<-3
plot(eModel, cex = .5, label.offset = .02, edge.width = ew)

#Bootstrap support
boot_sup<-l1ou_bootstrap_support(eModel, nItrs = 2, multicore = T, nCores = 10)
e.l<-round(boot_sup$detection.rate *100, digits = 1)
e.l<-ifelse(e.l>10, paste0(e.l, "%"), NA)
plot(eModel, edge.label = e.l, edge.ann.cex = .7, edge.label.ann = T, cex = .5,
             label.offset = .02, edge.width = ew)

saveRDS(boot_sup, "output/eModel_boot_sup.rds")

# Constrained set of candidate branches with a shift ####
eModel_c<-estimate_shift_configuration(lizard$tree, lizard$Y, criterion = "AICc", nCores = 3)
ce<-eModel_c$shift.configuration
eModel_cd<-estimate_shift_configuration(lizard$tree, lizard$Y, candid.edges = ce, nCores = 3)
plot(eMode_cd, edge.ann.cex = .7, cex = .5, label.offset = .02)

saveRDS(eModel_c, "output/eModel_c.rds")
saveRDS(eModel_cd, "output/eModel_cd.rds")

# Fit OU model given a configuration ####
eModel_ou<-fit_OU(eModel$tree, eModel$Y, eModel$profile$configurations[[2]],
                  l1ou.options = eModel$l1ou.options)
saveRDS(eModel_ou, "output/eModel_ou.rds")
plot(eModel_ou)
