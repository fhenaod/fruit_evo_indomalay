library(dplyr)
library(treeplyr)
library(bayou)

# load & prepare data ####
spptree1<-read.tree("data/zanne_tree.tre")
spptraits1<-read.csv("data/spptraits1.csv")

tree_data<-make.treedata(spptree1, spptraits1)
tree_data_ln<-mutate(tree_data, ln_fruit_lg = log(fruit_lg), ln_seed_lg = log(seed_lg), ln_num_seeds = log(num_seeds))
sd2<-sqrt(log(1+ (var(pull(tree_data_ln$dat,fruit_lg), na.rm = T)) / (mean(pull(tree_data_ln$dat,fruit_lg), na.rm = T))^2 ))

d_fruit_lg_ln<-(filter(tree_data_ln, !is.na(ln_fruit_lg)))
apply(select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku, newguinea), 2, function(x) sum(is.na(x)))
d_fruit_lg_ln <- filter(d_fruit_lg_ln, !is.na(sunda), !is.na(sulawesi), !is.na(maluku), !is.na(newguinea))

cache<-bayou:::.prepare.ou.univariate(d_fruit_lg_ln$phy, d_fruit_lg_ln[["ln_fruit_lg"]], pred = select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku, newguinea))
pred <- cache$pred
.prepare.ou.univariate <- bayou:::.prepare.ou.univariate
.tipregime <- bayou:::.tipregime

#cutoff <- .2
#sumpars <- list(sb = which(sum_c$branch.posteriors$pp > cutoff))
#sumpars$k <- length(sumpars$sb)
#sumpars$ntheta <- length(sumpars$sb)+1
#sumpars$loc <- rep(0, sumpars$k)
#sumpars$t2 <- 2:sumpars$ntheta

## priors ####
par.alpha <- list(scale = 1)
par.sig2 <- list(scale = 1)
par.beta_sunda     <- list(mean = 0, sd = .1) 
par.beta_sulawesi  <- list(mean = 0, sd = .1) 
par.beta_maluku    <- list(mean = 0, sd = .1) 
par.beta_newguinea <- list(mean = 0, sd = .1) 

par.k <- list(lambda = (2*Ntip(d_fruit_lg_ln$phy)-2)*(2.5/100), kmax = (2*Ntip(d_fruit_lg_ln$phy)-2)*(5/100))
par.sb <- list(bmax = 1, prob = 1)
par.theta <- list(mean = mean(getVector(d_fruit_lg_ln, ln_fruit_lg)), sd = sd2)

bmax <- prob <- rep(1, nrow(cache$edge)) 

# Tuning pars & model making ####
D.XXX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_sunda = 0.1, beta_sulawesi = 0.004, beta_maluku = 0.2, beta_newguinea = 0.05, k = rep(1, nrj), theta = 3,    slide = 1, missing.pred = 1)
D.1XX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_sunda = 0.1, beta_sulawesi = 0.004, beta_maluku = 0.2, beta_newguinea = 0.05, k = rep(1, nrj), theta = 0.15, slide = 1, missing.pred = 1)
D.XNX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_sunda = 0.2, beta_sulawesi = 0.004, beta_maluku = 0.2, beta_newguinea = 0.05, k = rep(1, nrj), theta = 1,    slide = 1, missing.pred = 1)
# .2 -.4 if out pars are more than .4 it should be move up, if lower tunning lower. They should match up.

# Code explanation. 5 numbers given either: R - RJ; N - Fixed multiple; 1 - Fixed global; 0 - Absent
# 1:β0; 2: β_sunda; 3: β_sulawesi; 4: β_maluku; 5: β_newguinea 

gens <- 10000

## RJ intercept & slope ####
# intercept and sunda
prior.RR000 <- make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                          dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                       dbeta_sunda = "dnorm",
                                       #dbeta_sulawesi = "dnorm",
                                       #dbeta_maluku = "dnorm",
                                       #dbeta_newguinea = "dnorm",
                                       dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                          param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                       dbeta_sunda = par.beta_sunda, sd = sd2, # regresion coeff. "dbeta_predictor"
                                       dk = par.k, dsb = par.sb, dtheta = par.theta)
                          )

model.RR000 <- makeBayouModel(ln_fruit_lg ~ sunda, rjpars = c("theta", "sunda"), tree = d_fruit_lg_ln$phy,  
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg),
                              pred = d_fruit_lg_ln$dat, prior.RR000, D = D.XXX(2))

prior.RR000(model.RR000$startpar)
model.RR000$model$lik.fn(model.RR000$startpar, cache, cache$dat)$loglik

mcmc.RR000 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda)), 
                       model = model.RR000, prior = prior.RR000, startpar = model.RR000$startpar,
                       new.dir = TRUE, outname = paste(mod, "r1", sep = "_"), plot.freq = NULL, 
                       ticker.freq = 2000, samp = 200)
mcmc.RR000$run(gens)

## reversable intercept fixed sunda
prior.R1000 <- make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                          dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                       dbeta_sunda = "dnorm",
                                       #dbeta_sulawesi = "dnorm",
                                       #dbeta_maluku = "dnorm",
                                       #dbeta_newguinea = "dnorm",
                                       dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                          param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                       dbeta_sunda = par.beta_sunda, sd = sd2, # regresion coeff. "dbeta_predictor"
                                       dk = par.k, dsb = par.sb, dtheta = par.theta)
)

model.R1000 <- makeBayouModel(ln_fruit_lg ~ sunda, rjpars = c("theta"), tree = d_fruit_lg_ln$phy,  
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg),
                              pred = d_fruit_lg_ln$dat, prior.R1000, D = D.XXX(1))

prior.R1000(model.R1000$startpar)
model.R1000$model$lik.fn(model.R1000$startpar, cache, cache$dat)$loglik

mcmc.R1000 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                             pred = select(d_fruit_lg_ln$dat, c(sunda)), 
                             model = model.R1000, prior = prior.R1000, startpar = model.R1000$startpar,
                             new.dir = TRUE, outname = paste(mod, "r1", sep = "_"), plot.freq = NULL, 
                             ticker.freq = 2000, samp = 200)
mcmc.R1000$run(gens)

## Simple regression ####
# sunda
prior.11000 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                        dbeta_sunda = "dnorm",
                                        #dbeta_sulawesi = "dnorm",
                                        #dbeta_maluku = "dnorm",
                                        #dbeta_newguinea = "dnorm",
                                        dsb = "fixed", dk = "fixed", dtheta = "dnorm"), 
                           param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                        dbeta_sunda = par.beta_sunda,
                                        #dbeta_sulawesi = par.beta_sulawesi,
                                        #dbeta_maluku = par.beta_maluku,
                                        #dbeta_newguinea = par.beta_newguinea,
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.11000 <- makeBayouModel(ln_fruit_lg ~ sunda, rjpars = numeric(0), tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda),
                              prior = prior.11000, impute = NULL, D = D.1XX(1))

prior.11000(model.11000$startpar)
model.11000$model$lik.fn(model.11000$startpar, cache, cache$dat)$loglik

mcmc.11000 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda)), 
                         model = model.11000$model, prior = prior.11000, startpar = model.11000$startpar, 
                         new.dir = "mod_11000/", outname = "11000_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
saveRDS(mcmc.11000, file = "mcmc.11000.rds")

mcmc.11000$run(gens)
chain.N1100 <- mcmc.11000$load()
chain.11000 <- set.burnin(chain.11000, .3)
saveRDS(chain.11000, file = "chain.11000.rds")
sum_c11000 <- summary(chain.11000)
saveRDS(sum_c11000, file = "sum.11000.rds")

require(foreach)
require(doParallel)
registerDoParallel(cores = 25)
Bk <- qbeta(seq(0,1, length.out = 30),.3, 1)
ss.11000 <- mcmc.11000$steppingstone(gens, chain.11000, Bk = Bk, plot = F)
saveRDS(ss.11000, file = "ss.11000.rds")
print(ss.11000$lnr)

# sunda + sulawesi 
prior.11100 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                        dbeta_sunda = "dnorm",
                                        dbeta_sulawesi = "dnorm",
                                        #dbeta_maluku = "dnorm",
                                        #dbeta_newguinea = "dnorm",
                                        dsb = "fixed", dk = "fixed", dtheta = "dnorm"), 
                           param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                        dbeta_sunda = par.beta_sunda,
                                        dbeta_sulawesi = par.beta_sulawesi,
                                        #dbeta_maluku = par.beta_maluku,
                                        #dbeta_newguinea = par.beta_newguinea,
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.11100 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                            dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda, sulawesi),
                            prior = prior.11100, impute = NULL, D = D.1XX(1))

prior.11100(model.11100$startpar)
model.11100$model$lik.fn(model.11100$startpar, cache, cache$dat)$loglik

mcmc.11100 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi)), 
                         model = model.11100$model, prior = prior.11100, startpar = model.11100$startpar, 
                         new.dir = "mod_11100/", outname = "11100_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mcmc.11100$run(gens)

# sunda + sulawesi + maluku 
prior.11110 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                        dbeta_sunda = "dnorm",
                                        dbeta_sulawesi = "dnorm",
                                        dbeta_maluku = "dnorm",
                                        #dbeta_newguinea = "dnorm",
                                        dsb = "fixed", dk = "fixed", dtheta = "dnorm"), 
                           param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                        dbeta_sunda = par.beta_sunda,
                                        dbeta_sulawesi = par.beta_sulawesi,
                                        dbeta_maluku = par.beta_maluku,
                                        #dbeta_newguinea = par.beta_newguinea,
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.11110 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi + maluku, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku),
                              prior = prior.11110, impute = NULL, D = D.1XX(1))

prior.11110(model.11110$startpar)
model.11110$model$lik.fn(model.11110$startpar, cache, cache$dat)$loglik

mcmc.11110 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi, maluku)), 
                         model = model.11110$model, prior = prior.11110, startpar = model.11110$startpar, 
                         new.dir = "mod_11110/", outname = "11110_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mcmc.11110$run(gens)

# sunda + sulawesi + maluku + newguinea
prior.11111 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                        dbeta_sunda = "dnorm",
                                        dbeta_sulawesi = "dnorm",
                                        dbeta_maluku = "dnorm",
                                        dbeta_newguinea = "dnorm",
                                        dsb = "fixed", dk = "fixed", dtheta = "dnorm"), 
                           param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                        dbeta_sunda = par.beta_sunda,
                                        dbeta_sulawesi = par.beta_sulawesi,
                                        dbeta_maluku = par.beta_maluku,
                                        dbeta_newguinea = par.beta_newguinea,
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.11111 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi + maluku + newguinea, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku, newguinea),
                              prior = prior.11111, impute = NULL, D = D.1XX(1))

prior.11111(model.11111$startpar)
model.11111$model$lik.fn(model.11111$startpar, cache, cache$dat)$loglik

mcmc.11111 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi, maluku, newguinea)), 
                         model = model.11111$model, prior = prior.11111, startpar = model.11111$startpar, 
                         new.dir = "mod_11111/", outname = "11111_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mcmc.11111$run(gens)

# newguinea 
prior.10001 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE,
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy",
                                        # dbeta_sunda = "dnorm",
                                        #dbeta_sulawesi = "dnorm",
                                        #dbeta_maluku = "dnorm",
                                        dbeta_newguinea = "dnorm",
                                        dsb = "fixed", dk = "fixed", dtheta = "dnorm"),
                           param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                        # dbeta_sunda = par.beta_sunda,
                                        # dbeta_sulawesi = par.beta_sulawesi,
                                        #dbeta_maluku = par.beta_maluku,
                                        dbeta_newguinea = par.beta_newguinea,
                                        dk = "fixed", dsb = "fixed",
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.10001 <- makeBayouModel(ln_fruit_lg ~ newguinea, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg),
                              pred = select(d_fruit_lg_ln$dat, newguinea),
                              prior = prior.10001, impute = NULL, D = D.1XX(1))

prior.10001(model.10001$startpar)
model.10001$model$lik.fn(model.10001$startpar, cache, cache$dat)$loglik

mcmc.10001 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg),
                             pred = select(d_fruit_lg_ln$dat, c( newguinea)),
                             model = model.10001$model, prior = prior.10001, startpar = model.10001$startpar,
                             new.dir = "mod_10001/", outname = "10001_r1", plot.freq = NULL,
                             ticker.freq = 2000, samp = 200, perform.checks = T)
mcmc.10001$run(gens)

## Separate slope ####
# sunda
prior.1N000 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                        dbeta_sunda = "dnorm",
                                        #dbeta_sulawesi = "dnorm",
                                        #dbeta_maluku = "dnorm",
                                        #dbeta_newguinea = "dnorm",
                                        dsb = "fixed", dk = "fixed", dtheta = "dnorm"), 
                           param = list(dalpha = par.alpha, dsig2 = par.sig2,
                                        dbeta_sunda = par.beta_sunda,
                                        #dbeta_sulawesi = par.beta_sulawesi,
                                        #dbeta_maluku = par.beta_maluku,
                                        #dbeta_newguinea = par.beta_newguinea,
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.1N000 <- makeBayouModel(ln_fruit_lg ~ sunda, rjpars = c("sunda"), tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg),
                              pred = select(d_fruit_lg_ln$dat, sunda),
                              prior = prior.1N000, impute = NULL, D = D.XNX(1))

prior.1N000(model.1N000$startpar)
model.1N000$model$lik.fn(model.1N000$startpar, cache, cache$dat)$loglik

mcmc.1N000 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda)), 
                         model = model.1N000$model, prior = prior.1N000, startpar = model.1N000$startpar, 
                         new.dir = "mod_1N000/", outname = "1N000_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mcmc.1N000$run(gens)
