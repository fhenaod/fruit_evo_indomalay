library(treeplyr)
library(bayou)
library(phytools)

## load & prepare data ####
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

sumpars <- readRDS("chain2/model_free/shiftsum_free.rds")
sum_c <- readRDS("chain2/chain_sum.rds")
#cutoff <- .3
#sumpars$sb <- list(sb = which(sum_c$branch.posteriors$pp > cutoff))
fixed.k <- sumpars$pars$k 
fixed.sb <- sumpars$pars$sb
fixed.loc <- rep(0, sumpars$pars$k)
sumpars$t2 <- 2:sumpars$pars$ntheta

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

## Tuning pars & model making ####
D.XXX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_sunda = 0.1, beta_sulawesi = 0.004, beta_maluku = 0.2, beta_newguinea = 0.05, k = rep(1, nrj), theta = 3,    slide = 1, missing.pred = 1)
D.1XX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_sunda = 0.1, beta_sulawesi = 0.004, beta_maluku = 0.2, beta_newguinea = 0.05, k = rep(1, nrj), theta = 0.15, slide = 1, missing.pred = 1)
D.XNX <- function(nrj) list(alpha = 0.7, sig2 = 0.5, beta_sunda = 0.2, beta_sulawesi = 0.004, beta_maluku = 0.2, beta_newguinea = 0.05, k = rep(1, nrj), theta = 1,    slide = 1, missing.pred = 1)
# .2 -.4 if out pars are more than .4 it should be move up, if lower tunning lower. They should match up.

# Code explanation. 5 numbers given either: R - RJ; N - Fixed multiple; 1 - Fixed global; 0 - Absent
# 1:β0; 2: β_sunda; 3: β_sulawesi; 4: β_maluku; 5: β_newguinea 

gens <- 10000

## RJ intercept & slope ####
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
                             pred = select(d_fruit_lg_ln$dat, c(newguinea)),
                             model = model.10001$model, prior = prior.10001, startpar = model.10001$startpar,
                             new.dir = "mod_10001/", outname = "10001_r1", plot.freq = NULL,
                             ticker.freq = 2000, samp = 200, perform.checks = T)
mcmc.10001$run(gens)

## Separate  ####
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

# separate intercepts
#sunda
prior.N1000 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = F, 
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
                           fixed = list(k = fixed.k, sb = fixed.sb, loc = fixed.loc)
)

model.N1000 <- makeBayouModel(ln_fruit_lg ~ sunda, rjpars = c("theta"), tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg),
                              pred = select(d_fruit_lg_ln$dat, sunda),
                              prior = prior.N1000, impute = NULL, D = D.XXX(1))

prior.N1000(model.N1000$startpar)
model.N1000$model$lik.fn(model.N1000$startpar, cache, cache$dat)$loglik

mcmc.N1000 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                             pred = select(d_fruit_lg_ln$dat, c(sunda)), 
                             model = model.N1000$model, prior = prior.N1000, startpar = model.N1000$startpar, 
                             new.dir = "mod_N1000/", outname = "N1000_r1", plot.freq = NULL, 
                             ticker.freq = 2000, samp = 200, perform.checks = T)

#sunda + sulawesi + maluku + new_guinea
prior.N1111 <-  make.prior(d_fruit_lg_ln$phy, plot.prior = F, 
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
                           fixed = list(k = fixed.k, sb = fixed.sb, loc = fixed.loc)
)

model.N1111 <- makeBayouModel(ln_fruit_lg ~ sunda + sulawesi + maluku + newguinea, rjpars = c("theta"), tree = d_fruit_lg_ln$phy,
                              dat = getVector(d_fruit_lg_ln, ln_fruit_lg),
                              pred = select(d_fruit_lg_ln$dat, sunda, sulawesi, maluku, newguinea),
                              prior = prior.N1111, impute = NULL, D = D.XXX(1))

prior.N1111(model.N1111$startpar)
model.N1111$model$lik.fn(model.N1111$startpar, cache, cache$dat)$loglik

mcmc.N1111 <- bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, ln_fruit_lg), 
                             pred = select(d_fruit_lg_ln$dat, c(sunda, sulawesi, maluku, newguinea)), 
                             model = model.N1111$model, prior = prior.N1111, startpar = model.N1111$startpar, 
                             new.dir = "mod_N1111/", outname = "N1111_r1", plot.freq = NULL, 
                             ticker.freq = 2000, samp = 200, perform.checks = T)

## Stochastic maping ####
sunda_ard <- make.simmap(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, sunda), model = "ARD", nsim = 1000)
sula_ard <- make.simmap(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, sulawesi), model = "ARD", nsim = 1000)
malu_ard <- make.simmap(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, maluku), model = "ARD", nsim = 1000)
new_ard <- make.simmap(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, newguinea), model = "ARD", nsim = 1000)

sunda_ard <- readRDS("simmap/sunda_ard.rds")
sula_ard <- readRDS("simmap/sula_ard.rds")
malu_ard <- readRDS("simmap/malu_ard.rds")
new_ard <- readRDS("simmap/new_ard.rds")

sunda_sum <- summary(sunda_ard, plot = F)
sula_sum <- summary(sula_ard, plot = F)
malu_sum <- summary(malu_ard, plot = F)
new_sum <- summary(new_ard, plot = F)

cols <- setNames(palette()[1:length(levels(factor(getVector(tree_data, maluku))))], levels(factor(getVector(tree_data, maluku))))
plot(new_sum, cols, fsize = .01, lwd = .5, cex = c(.4,.2))
add.simmap.legend(colors = cols, x = 0.95*par()$usr[1], y = 0.9*par()$usr[4], prompt = FALSE, fsize = 0.8)

sunda_dens_m <- densityMap(sula_ard)

sunda_dens_m <- readRDS("simmap/sunda_densM.rds")
sula_dens_m <- readRDS("simmap/sula_densM.rds")
malu_dens_m <- readRDS("simmap/malu_densM.rds")
new_dens_m <- readRDS("simmap/new_densM.rds")

plot(sunda_dens_m, fsize = c(0.1, 1), lwd = .5)

disp_er_sum <- summary(disp_simp_er, plot = F)
disp_ard_sum <- summary(disp_simp_ard, plot = F)
disp_sym_sum <- summary(disp_simp_sym, plot = F)

saveRDS(disp_er_sum , "output/disp_er_sum.rds")
saveRDS(disp_ard_sum, "output/disp_ard_sum.rds")
saveRDS(disp_sym_sum, "output/disp_sym_sum.rds")
