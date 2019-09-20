library(dplyr)
library(treeplyr)
library(bayou)

# Load data ####
spptree1<-read.tree("data/spptreeZanneShort.tree") 
spptraits1<-read.csv("data/spptraits1.csv")

par(mfrow = c(3,1))
hist(log(spptraits1$fruit_lg), breaks = 20, main = "", xlab = "Log Fruit length", ylim = c(0, 500), xlim = c(0,7))
abline(v = log(mean(spptraits1$fruit_lg, na.rm = T)), col = "red", lwd = 2, lty = 2)

hist(log(spptraits1$seed_lg), breaks = 20, main = "", xlab = "Log Seed length", ylim = c(0, 150), xlim = c(-1,6))
abline(v = log(mean(spptraits1$seed_lg, na.rm = T)), col = "red", lwd = 2, lty = 2)

hist(log(spptraits1$num_seeds), breaks = 20, main = "", xlab = "Log Number of Seeds", ylim = c(0, 250), xlim = c(0,6))
abline(v = log(mean(spptraits1$num_seeds, na.rm = T)), col = "red", lwd = 2, lty = 2)

# make tree_data
tree_data<-make.treedata(spptree1, spptraits1)
tree_data_ln<-mutate(tree_data, ln_fruit_lg = log(fruit_lg), ln_seed_lg = log(seed_lg), ln_num_seeds = log(num_seeds))

d_fruit_lg_ln_ln<-(filter(tree_data_ln, !is.na(ln_fruit_lg)))
name.check(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln_ln, ln_fruit_lg))

# Prior & rjMCMC run ####

# OU free/global model
priorOU<-make.prior(d_fruit_lg_ln$phy, 
                      dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                   dk = "cdpois", dtheta = "dnorm"),
                      param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                   dk = list(lambda = 10, kmax = Ntip(d_fruit_lg_ln$phy)-1), dsb = list(bmax = 1, prob = 1), 
                                   dtheta = list(mean = mean(getVector(d_fruit_lg_ln, fruit_lg)), sd = 1.5*sd(getVector(d_fruit_lg_ln, fruit_lg))))
)

quantiles<-c(0, 0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99, 1)
alfs<-rhalfcauchy(10000, scale = 0.1)
half_l<-log(2/alfs)
qs<-quantile(half_l, quantiles) 
round(qs, 2)

startpars<-priorSim(priorOU, d_fruit_lg_ln$phy, plot=TRUE)$pars[[1]]
priorOU(startpars)

mcmcOU<-bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, fruit_lg),
                       prior = priorOU, new.dir = "modelOU/", outname = "modelOU_r001", plot.freq = NULL) 
mcmcOU$run(10000)

chainOU<-mcmcOU$load(saveRDS = T, file = 'modelOU/mcmcOU.RDS')
chainOU<-set.burnin(chainOU, 0.3)
summary(chainOU)
plot(chainOU, auto.layout = FALSE)

plotSimmap.mcmc(chainOU, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(d_fruit_lg_ln$phy, chainOU, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, fruit_lg), burnin = 0.3, chainOU, pp.cutoff = 0.3)

## Blunderburst model (BB)
par.halflife<-list(meanlog=2.5, sdlog=2.5)
par.Vy<-list(meanlog=log(0.0958), sdlog=0.25)

priorBB<-make.prior(d_fruit_lg_ln$phy, 
                      dists = list(dhalflife = "dlnorm", dVy = "dlnorm", 
                                   dk = "cdpois", dsb = "dsb", dtheta = "dnorm"),
                      param = list(dhalflife = par.halflife,
                                   dVy = par.Vy,
                                   dk = list(lambda = 10, kmax = 50), dsb = list(bmax = 1, prob = 1), 
                                   dtheta = list(mean = getVector(d_fruit_lg_ln, fruit_lg), sd = 1.5*sd(getVector(d_fruit_lg_ln, fruit_lg)))),
                      model = "OUrepar"
)

mcmcBB<-bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, fruit_lg), model = "OUrepar", 
                       prior = priorBB, new.dir = "modelBB/", outname = "modelBB_r001", plot.freq = NULL)
mcmcBB$run(10000)
chainBB<-mcmcBB$load(saveRDS = T, file = 'modelBB/mcmcBB.RDS')
chainBB<-set.burnin(chainBB, 0.3)
summary(chainBB)
plot(chainBB)

# Quantitative Genetics model (QG)
par.h2<-list(shape1=10, shape2=10)
par.P<-list(meanlog=log(0.12), sdlog=0.2)
par.w2<-list(meanlog=log(100), sdlog=2.5)
par.Ne<-list(meanlog=log(500000), sdlog=2.5)

QGtree<-d_fruit_lg_ln$phy
QGtree$edge.length<-QGtree$edge.length/2 ## Convertion to generation times

priorQG<-make.prior(QGtree, 
                      dists = list(dh2 = "dbeta", dP = "dlnorm",
                                   dw2 = "dlnorm", dNe = "dlnorm",
                                   dk = "cdpois", dtheta = "dnorm"),
                      param = list(dh2 = par.h2,
                                   dP = par.P,
                                   dw2 = par.w2,
                                   dNe = par.Ne,
                                   dk = list(lambda = 10, kmax = 50), dsb = list(bmax = 1, prob = 1), 
                                   dtheta = list(mean = mean(getVector(d_fruit_lg_ln, fruit_lg)), sd = 1.5*sd(getVector(d_fruit_lg_ln, fruit_lg)))),
                      model = "QG"
)

mcmcQG<-bayou.makeMCMC(QGtree, getVector(d_fruit_lg_ln, fruit_lg), model = "QG", 
                       startpar = NULL, prior = priorQG, new.dir="modelQG/", outname="modelQG_r001", plot.freq = NULL)
mcmcQG$run(10000)

chainQG<-mcmcQG$load(saveRDS = T, file = 'modelQG/mcmcQG.RDS')
chainQG<-set.burnin(chainQG, 0.3)
summary(chainQG)
plot(chainQG)

par(mfrow=c(2,2))
plotBayoupars(truepars, tree, main = "True parameters")
plotSimmap.mcmc(chainQG, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(tree, chainQG, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(tree, dat, burnin = 0.3, chainQG, pp.cutoff = 0.3)

# Model comparison ####
registerDoParallel(cores = 3)
Bk<-qbeta(seq(0,1, length.out=5), 0.3,1)
ssOU<-mcmcOU$steppingstone(10000, chainOU, Bk, burnin = 0.3, plot = FALSE)
ssBB<-mcmcBB$steppingstone(10000, chainBB, Bk, burnin = 0.3, plot = FALSE)
ssQG<-mcmcQG$steppingstone(10000, chainQG, Bk, burnin = 0.3, plot = FALSE)
mlnL<-c("OU" = ssOU$lnr, "BB" = ssBB$lnr, "QG" = ssQG$lnr)
mlnL

plot(ssOU)
plot(ssBB)
plot(ssQG)

# Fixed models ####

# Fixed prior, an hypothesis
trueFixedPrior<-make.prior(d_fruit_lg_ln$phy, 
                           dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                        dk = "fixed", dsb = "fixed", dtheta = "dnorm"),
                           param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                        dk = "fixed", dsb = "fixed", 
                                        dtheta = list(mean = mean(getVector(d_fruit_lg_ln, fruit_lg)), sd = .5*sd(getVector(d_fruit_lg_ln, fruit_lg)))),
                           fixed = list(k = truepars$k, sb = truepars$sb)
)

# Fixed prior, alternative hypothesis
altlocations<-list(sb = c(89, 70, 85, 47, 50), loc = c(1.5, 6.6, 5.2, 3.7, 11.1))
altpars<-truepars ####
altpars$k<-5
altpars$sb<-altlocations$sb
altpars$loc<-altlocations$loc
altpars$t2<-c(2, 3, 4, 5, 3) # a covergence in regimes '3'

alternativeFixedPrior<-make.prior(d_fruit_lg_ln$phy, 
                                  dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", 
                                               dk = "fixed", dsb = "fixed", dtheta = "dnorm"),
                                  param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                               dk = "fixed", dsb = "fixed", 
                                               dtheta = list(mean = mean(getVector(d_fruit_lg_ln, fruit_lg)), sd = 1.5*sd(getVector(d_fruit_lg_ln, fruit_lg)))),
                                  fixed = list(k = 5, ntheta = 5, sb = altpars$sb, loc = altpars$loc, t2 = altpars$t2)
)

par(mfrow=c(1,2))
plotBayoupars(truepars, tree, main="True Pars")
plotBayoupars(altpars, tree, main="Alternative Hypothesis")

mcmcFixed1<-bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, fruit_lg), SE = MEvar, prior = trueFixedPrior, 
                          new.dir = "Fixed/", outname = "modelTrueFixed_r001", plot.freq = NULL)
mcmcFixed1$run(10000)

mcmcFixed2<-bayou.makeMCMC(d_fruit_lg_ln$phy, getVector(d_fruit_lg_ln, fruit_lg), SE = MEvar, prior = alternativeFixedPrior, 
                           new.dir = "Fixed/", outname = "modelAltFixed_r001", plot.freq = NULL)
mcmcFixed2$run(10000)

chainFixed1<-mcmcFixed1$load(saveRDS = T, file = 'Fixed/mcmcFixed1.RDS')
chainFixed2<-mcmcFixed2$load(saveRDS = T, file = 'Fixed/mcmcFixed2.RDS')

# Marginal Likelihood estimation
Bk<-qbeta(seq(0,1, length.out=5), 0.3,1)
ssFixed1<-mcmcFixed1$steppingstone(10000, chainFixed1, Bk)
ssFixed2<-mcmcFixed2$steppingstone(10000, chainFixed2, Bk)
ssFixed1$lnr
ssFixed2$lnr

# Customized models ####

# Global intercepts & slopes
prior.11<-make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                     dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dbeta_lnMass = "dnorm",
                                  dsb = "fixed", dk = "fixed", dtheta = "dnorm"), 
                     param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                  dbeta_lnMass = list(mean = 0.7, sd = 0.15), # regresion coeff. "dbeta_predictor"
                                  dtheta = list(mean = 0, sd = 1)),
                     fixed = list(k = 0, sb = numeric(0))
)

# Separate intercepts & global slope
prior.N1<-make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                     dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dbeta_lnMass = "dnorm",
                                  dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                     param = list(dbeta_lnMass = list(mean = 0.7, sd = 0.15), ## Why there are no alpha and dsig pars? ##
                                  dk = list(lambda = 10, kmax = 50),
                                  dtheta = list(mean = 0, sd = 1))
)

# Separate intercepts & slopes
prior.NN<-make.prior(d_fruit_lg_ln$phy, plot.prior = FALSE, 
                     dists = list(dalpha = "dhalfcauchy", dsig2 = "dhalfcauchy", dbeta_lnMass = "dnorm",
                                  dsb = "dsb", dk = "cdpois", dtheta = "dnorm"), 
                     param = list(dalpha = list(scale = 0.1), dsig2 = list(scale = 0.1),
                                  dbeta_lnMass = list(mean = 0.7, sd = 0.15),
                                  dk = list(lambda = 10, kmax = 50), 
                                  dtheta = list(mean = 0, sd = 1))
)

# Tuning pars & make models
D11<-list(alpha = 2, sig2 = 2, beta_lnMass = 0.1, k = 1,      theta = 0.5, slide = 1)
DN1<-list(alpha = 2, sig2 = 2, beta_lnMass = 0.1, k = 1,      theta = 2,   slide = 1)
DNN<-list(alpha = 2, sig2 = 2, beta_lnMass = 0.3, k = c(1,1), theta = 2,   slide = 1)

model.11<-makeBayouModel(dat ~ lnMass, rjpars = c(), 
                         tree = d_fruit_lg_ln$phy, dat = dat, pred = pred, SE = MEvar, prior = prior.11, D = D11)

model.N1<-makeBayouModel(dat ~ lnMass, rjpars = c("theta"),  
                         tree = d_fruit_lg_ln$phy, dat = dat, pred = pred, SE = MEvar, prior = prior.N1, D = DN1)

model.NN<-makeBayouModel(dat ~ lnMass, rjpars = c("theta", "lnMass"),  
                         tree = d_fruit_lg_ln$phy, dat = dat, pred = pred, SE = MEvar, prior = prior.NN, D = DNN)

mcmc.11<-bayou.makeMCMC(d_fruit_lg_ln$phy, dat, pred=pred, SE=MEvar, model=model.11$model, prior=prior.11, 
                          startpar=model.11$startpar, new.dir="Allometry/", outname="model11_r001", plot.freq=NULL)
mcmc.N1<-bayou.makeMCMC(d_fruit_lg_ln$phy, dat, pred=pred, SE=MEvar, model=model.N1$model, prior=prior.N1, 
                          startpar=model.N1$startpar, new.dir="Allometry/", outname="modelN1_r001", plot.freq=NULL)
mcmc.NN<-bayou.makeMCMC(d_fruit_lg_ln$phy, dat, pred=pred, SE=MEvar, model=model.NN$model, prior=prior.NN, 
                          startpar=model.NN$startpar, new.dir="Allometry/", outname="modelNN_r001", plot.freq=NULL)

mcmc.11$run(10000)
mcmc.N1$run(10000)
mcmc.NN$run(10000)

chain.11<-set.burnin(mcmc.11$load(saveRDS = T, file = 'Allometry/mcmc11.RDS'), 0.3)
chain.N1<-set.burnin(mcmc.N1$load(saveRDS = T, file = 'Allometry/mcmcN1.RDS'), 0.3)
chain.NN<-set.burnin(mcmc.NN$load(saveRDS = T, file = 'Allometry/mcmcNN.RDS'), 0.3)

shiftsumsN1<-shiftSummaries(chain.N1, mcmc.N1, pp.cutoff = 0.5)
shiftsumsNN<-shiftSummaries(chain.NN, mcmc.NN, pp.cutoff = 0.5)

plotShiftSummaries(shiftsumsN1, lwd = 2, single.plot = TRUE, label.pts = FALSE)
plotShiftSummaries(shiftsumsNN, lwd = 2, single.plot = TRUE, label.pts = FALSE)

registerDoParallel(cores = 3)
Bk<-qbeta(seq(0,1, length.out = 5), 0.3,1)
ss.11<-mcmc.11$steppingstone(10000, chain.11, Bk, burnin = 0.3, plot = FALSE)
ss.N1<-mcmc.N1$steppingstone(10000, chain.N1, Bk, burnin = 0.3, plot = FALSE)
ss.NN<-mcmc.NN$steppingstone(10000, chain.NN, Bk, burnin = 0.3, plot = FALSE)
mlnL<-c("11" = ss.11$lnr, "N1" = ss.N1$lnr, "NN" = ss.NN$lnr)
mlnL
