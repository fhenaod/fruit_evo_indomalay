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

## rj inter & sunda: RR000 ####
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

mod <- "RR000"
gens <- 10000
mymcmc <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), pred = d_fruit_lg_ln$dat, 
                       model = model.RR000, prior = prior.RR000, startpar = model.RR000$startpar,
                       new.dir = TRUE, outname = paste(mod, "run1", sep = "_"), plot.freq = NULL, 
                       ticker.freq = 2000, samp = 200)
mymcmc$run(gens)

## Separate intercepts ####
gens <- 10000
# sunda
prior.N1000 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
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
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.N1000 <- makeBayouModel(ln_fruit_lg ~ sunda, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda),
                              prior = prior.N1000, impute = NULL, D = D.XXX(1))

prior.N1000(model.N1000$startpar)
model.N1000$model$lik.fn(model.N1000$startpar, cache, cache$dat)$loglik

mymcmc <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda)), 
                         model = model.N1000$model, prior = prior.N1000, startpar = model.N1000$startpar, 
                         new.dir = "mod_N1000/", outname = "N1000_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mymcmc$run(gens)


# sunda + sulawesi 
prior.N1100 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
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
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = par.theta),
                           fixed = list(k = 0, sb = numeric(0), t2=numeric(0), loc = numeric(0))
)

model.N1100 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                            dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda, sulawesi),
                            prior = prior.N1100, impute = NULL, D = D.XXX(1))

prior.N1100(model.N1100$startpar)
model.N1100$model$lik.fn(model.N1100$startpar, cache, cache$dat)$loglik

mymcmc <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi)), 
                         model = model.N1100$model, prior = prior.N1100, startpar = model.N1100$startpar, 
                         new.dir = "mod_N1100/", outname = "N1100_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mymcmc$run(gens)

# sunda + sulawesi + maluku 
prior.N1110 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
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

model.N1110 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi + maluku, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku),
                              prior = prior.N1110, impute = NULL, D = D.XXX(1))

prior.N1110(model.N1110$startpar)
model.N1110$model$lik.fn(model.N1110$startpar, cache, cache$dat)$loglik

mymcmc <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi, maluku)), 
                         model = model.N1110$model, prior = prior.N1110, startpar = model.N1110$startpar, 
                         new.dir = "mod_N1110/", outname = "N1110_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mymcmc$run(gens)
# sunda + sulawesi + maluku + newguinea
prior.N1111 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
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

model.N1111 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi + maluku + newguinea, rjpars = NULL, tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg), pred = select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku, newguinea),
                              prior = prior.N1111, impute = NULL, D = D.XXX(1))

prior.N1111(model.N1111$startpar)
model.N1111$model$lik.fn(model.N1111$startpar, cache, cache$dat)$loglik

mymcmc <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                         pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi, maluku, newguinea)), 
                         model = model.N1111$model, prior = prior.N1111, startpar = model.N1111$startpar, 
                         new.dir = "mod_N1111/", outname = "N1111_r1", plot.freq = NULL, 
                         ticker.freq = 2000, samp = 200, perform.checks = T)
mymcmc$run(gens)