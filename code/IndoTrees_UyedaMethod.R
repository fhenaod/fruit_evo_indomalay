# Using the Uyuda et al. (2017; AmNat) method of modeling discrete transitions in regression coefficient across phylogenies
# Code modifed from https://github.com/ssb2017/bayou/blob/master/R/bayouWorkshopNotebook.Rmd


# Desktop
#setwd('D:/Box Sync/Projects/Projects (active)/Indo archipelago survey/Flora Malesiana analysis/MS Tree fruit evolution/Analysis')
# Laptop
#setwd('C:/Users/JB/Box Sync/Projects/Projects (active)/Indo archipelago survey/Flora Malesiana analysis/MS Tree fruit evolution/Analysis')


library(ape)
library(geiger)
library(phytools)
library(bayou)


#--------- DATA AND FUNCTIONS -----------------------------
rm(list = ls())
spptree1 <- ape::read.tree("spptreeZanneShort.tree") 
spptree1 <- reorder(spptree1, "postorder") # re-ordering the species in the tree
spptraits1 <- read.csv("spptraits1.csv")



# THIS IS UYEDA'S FUNCTION TO LOAD THE TRAIT DATA. I'M NOT SURE HOW TO MODIFY IT TO LOAD MY OWN DATA
# Now using the function *dataSim*, we can simulate trait data. 
dat <- dataSim(truepars, tree, model="OU")$dat
# To add realism, let's add some measurement error to the data. This is a good reminder to *always try to use measurement error in your analyses*. OU models especially are affected by measurement error. 
# This is because OU models have the effect of "erasing" evolutionary history with increasing *alpha*. If you don't account for measurement error, then that measurement error will be transferred to the 
# evolutionary process model. You can make a Brownian Motion model look very OU like if there is a lot of measurement error.
MEvar <- 0.1
dat <- dat + rnorm(length(dat), 0, sqrt(MEvar))



# Parameterization where priors are placed directly on phylogenetic half-life (*halflife*) and stationary variance (*Vy*), rather than *alpha* and *sig2*. 
# Here using a mildly informative prior on the phylogenetic half-life - a log-normal distribution:
par.halflife <- list(meanlog=2.5, sdlog=2.5)
par.Vy <- list(meanlog=log(0.0958), sdlog=0.25)


# Visualizing the priors
quantiles <- c(0, 0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99, 1)
samp <- rlnorm(10000, par.halflife$meanlog, par.halflife$sdlog)
hist(log(samp,10), breaks=100, main="Prior density of halflife")
abline(v=log(c(1,max(branching.times(tree))),10), col="red", lwd=2, lty=2)
# Notice that there is about equal density of prior probability on the half-life being greater than tree 
# height (rightmost red line) as there is below 1 million years (leftmost red line). The exact quantiles of this 
# distribution are:
qs <- qlnorm(quantiles, meanlog=par.halflife$meanlog, sdlog=par.halflife$sdlog)
round(setNames(qs, quantiles), 2)


# Make the prior 
priorBB <- make.prior(tree, dists=list(dhalflife="dlnorm", dVy="dlnorm", dk="cdpois", dsb="dsb", dtheta="dnorm"), param=list(dhalflife=par.halflife, dVy=par.Vy, dk=list(lambda=10, kmax=50), 
	dsb=list(bmax=1, prob=1), dtheta=list(mean=mean(dat), sd=1.5*sd(dat))), model="OUrepar")


# Make the MCMC object and run the chain
set.seed(1)
mcmcBB <- bayou.makeMCMC(tree, dat, SE=MEvar, model="OUrepar", prior=priorBB, new.dir="../output/", outname="modelBB_r001", plot.freq=NULL)
mcmcBB$run(10000)
chainBB <- mcmcBB$load()
chainBB <- set.burnin(chainBB, 0.3)
summary(chainBB)
plot(chainBB)


# Visualize the chain
par(mfrow=c(2,2))
plotBayoupars(truepars, tree, main = "True parameters")
plotSimmap.mcmc(chainBB, burnin = 0.3, pp.cutoff = 0.3)
plotBranchHeatMap(tree, chainBB, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(tree, dat, burnin = 0.3, chainBB, pp.cutoff = 0.3)
# Likely, we will have more shifts because we made the prior on *Vy* so narrow. Let's compare the posteriors from the two models.
quantile(chainOU$sig2/(2*chainOU$alpha), quantiles)
quantile(chainBB$Vy, quantiles)


# Model Comparison
# Alternative parameterizations, shift locations, and priors can be compared using Bayes Factors. This requires estimation of the marginal likelihood, which can be difficult. 
# **bayou** uses stepping-stone sampling to estimate the marginal likelihoods. To estimate marginal likelihoods, using the '$steppingstone' function in the mcmc object. For this exercise, 
# we will do a much shorter run than is recommended. If you have multiple cores available on your machine, you can make use of these to run the stepping stone analysis in parallel and 
# conduct the analysis much faster. 
Bk <- qbeta(seq(0,1, length.out=5), 0.3, 1)

# Marginal likelihood (increase the number of runs substantially for publication)
ssBB <- mcmcBB$steppingstone(10000, chainBB, Bk, burnin=0.3, plot=FALSE) # all sorts of info
mlBB <- ssBB$lnr # marginal likelihood

# If you get a couple errors it's probably OK, the algorithm takes the posterior and tries to fit various distributions to the parameters, and if it fails to optimize them it will throw an 
# error or two. However, as long as one of them fits OK it will run. Again, we have not run these for long enough or for enough steps (we prefer more like 50!), but you get the idea for how 
# you would proceed. Obviously, having more cores makes this go a LOT faster, and this is a computationally intensive procedure!
plot(ssBB)




#################################################################################

# 		Customized/Allometric models

# What if there is a known (or unknown) relationship between the trait of interest and another predictor variable? For example, we may be interested in a relationship between 
# a trait known to vary with body size, but consider the possibility that the relationship with body size itself varies over macroevolutionary time. Here, instead of having a
# single optimum that changes upon a regime shift, it is possible to have both the slope and intercept of the relationship change at once. bayou v2.0 allows you to include these
# additional predictors and test for shifts in the scaling between a trait and its predictors. 

# Let's simulate a dataset where the slope and intercept shift at different points in the tree. We're going to use the same shift locations as before, but add in a covariate
# with body size that changes in different parts of the tree. We also need to simulate the predictor data, in this case, let's use Brownian Motion.
set.seed(1)
tree <- sim.bdtree(b=1, d=0, stop="taxa", n=50, seed=1)
tree$edge.length <- tree$edge.length/max(branching.times(tree))*100
tree <- reorder(tree, "postorder")
truepars <- list(alpha = 0.5, sig2 = 0.05, k = 3, ntheta = 4, beta_lnMass = c(0.75, 0.6, 0.9, 0.67), theta = c(-1, 1.25, 0.5, 0), sb = c(94, 71, 50), loc = c(0, 0, 0), t2 = 2:4)
pred <- cbind("lnMass" = sim.char(tree, 0.2, model="BM", root=3)[,,1])
phytools::phenogram(tree, setNames(pred[,1], tree$tip.label), spread.labels=FALSE, main="Predictor: Body Size (lnMass)")
dat <- dataSim(truepars, tree, model="OU")$dat + truepars$beta_lnMass[bayou:::.tipregime(truepars, tree)] * pred[,1]


#------------------ JB NOTE: I don't understand how you input your own, real data into this, instead of simulating data. I can't even tell what format their "dat" object is in



# But our old visualization doesn't give the whole picture, because the trait covaries with body size:
par(mfrow=c(1,2))

# Plot the regime locations
plotRegimes(pars2simmap(truepars,tree)$tr,  col=pars2simmap(truepars,tree)$col)

# Plot the allometry
plot(pred[,1], dat, pch=21, bg=bayou:::.tipregime(truepars, tree), xlab="lnMass", "ylab"="Trait")

# Add the regression lines
dum <- lapply(1:truepars$ntheta, function(x) abline(truepars$theta[x], truepars$beta_lnMass[x],  lty=2, col=x))


# We are going to test 3 models in this analysis: Global intercepts & slopes (11), Separate intercepts & global slope (N1), and separate intercepts & slopes (NN).
# However, we're going to have to build these models to run them (they aren't built into **bayou**). As a convention, we're going to name our regression coefficients
# (other than the familiar intercept, *theta*) "beta_" followed by the predictor name (e.g. *beta_lnMass*). Here we imagine we have some fairly informative prior belief
# about what the allometry with body mass should be (normal distribution around 0.7).
prior.11 <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm", dsb="fixed", dk="fixed", dtheta="dnorm"), 
	param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1), dbeta_lnMass=list(mean=0.7, sd=0.15), dtheta=list(mean=0, sd=1)), fixed=list(k=0, sb=numeric(0)))

prior.N1 <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm", dsb="dsb", dk="cdpois", dtheta="dnorm"), 
	param=list(dbeta_lnMass=list(mean=0.7, sd=0.15), dk=list(lambda=10, kmax=50), dtheta=list(mean=0, sd=1)))

prior.NN <- make.prior(tree, plot.prior = FALSE, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_lnMass="dnorm", dsb="dsb", dk="cdpois", dtheta="dnorm"), 
	param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1), dbeta_lnMass=list(mean=0.7, sd=0.15), dk=list(lambda=10, kmax=50), dtheta=list(mean=0, sd=1)))

# Manually set tuning parameters, and make the models. There is a bit of art to tuning the parameters, which may require making multiple runs and trying to get the acceptance
# ratios in the right region (0.2-0.4). But these should work well for these models and data. If the acceptance ratio for a certain parameter is too high, increase the tuning
# parameter for that variable. If the acceptance ratio is too low, decrease it. The scale of the regression coefficient, for example, should give you some idea of what these
# parameters should be. 
D11 = list(alpha=2, sig2=2, beta_lnMass=0.1, k=1, theta=0.5, slide=1)
DN1 = list(alpha=2, sig2=2, beta_lnMass=0.1, k=1, theta=2, slide=1)
DNN = list(alpha=2, sig2=2, beta_lnMass=0.3, k=c(1,1), theta=2, slide=1)



# Now we use the function *makeBayouModel* to create a **bayou** model object that specifies all the components **bayou** needs to drop a new model into the analysis.
# Note that if you are interested in developing in **bayou** these are intended to be easy to make for customized models and there is a lot more possible than what is shown here.
# Note that each model only differs in the number of reversible-jump parameters (0, 1, and 2), the prior and the tuning parameters. By default, the starting regime map is plotted.

set.seed(1)
model.11 <- makeBayouModel(dat ~ lnMass, rjpars = c(), tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.11, D=D11)
model.N1 <- makeBayouModel(dat ~ lnMass, rjpars = c("theta"), tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.N1, D=DN1)
model.NN <- makeBayouModel(dat ~ lnMass, rjpars = c("theta", "lnMass"), tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.NN, D=DNN)

# We can now drop these model object into the analysis, along with the generated starting values (replacing the out of the box options of "OU", "OUrepar" and "QG"). 


# Make MCMC objects:
mcmc.11 <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.11$model, prior=prior.11, startpar=model.11$startpar, new.dir="../output/Allometry/", outname="model11_r001", plot.freq=NULL)
mcmc.N1 <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.N1$model, prior=prior.N1, startpar=model.N1$startpar, new.dir="../output/Allometry/", outname="modelN1_r001", plot.freq=NULL)
mcmc.NN <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.NN$model, prior=prior.NN, startpar=model.NN$startpar, new.dir="../output/Allometry/", outname="modelNN_r001", plot.freq=NULL)


# Run the models and load them in.
set.seed(1)
mcmc.11$run(10000)
mcmc.N1$run(10000)
mcmc.NN$run(10000)
chain.11 <- set.burnin(mcmc.11$load(), 0.3)
chain.N1 <- set.burnin(mcmc.N1$load(), 0.3)
chain.NN <- set.burnin(mcmc.NN$load(), 0.3)


# A particularly useful way to plot these is to use the *shiftSummaries* and *plotShiftSummaries* functions. Like other plotting functions, we define a posterior probability cutoff and only
# plot those shifts (*pp.cutoff*). Note that the global allometry (*11*), has no shifts and is not plotted here. 
shiftsumsN1 <- shiftSummaries(chain.N1, mcmc.N1, pp.cutoff=0.5, burnin=0.3)
shiftsumsNN <- shiftSummaries(chain.NN, mcmc.NN, pp.cutoff=0.5, burnin=0.3)
plotShiftSummaries(shiftsumsN1, lwd=2, single.plot=TRUE, label.pts=FALSE)
plotShiftSummaries(shiftsumsNN, lwd=2, single.plot=TRUE, label.pts=FALSE)


# As before, we can compare different models by estimating marginal likelihoods. Divide and conquer. 
registerDoParallel(cores=5)
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
set.seed(1)
ss.11 <- mcmc.11$steppingstone(10000, chain.11, Bk, burnin=0.3, plot=FALSE)
ss.N1 <- mcmc.N1$steppingstone(10000, chain.N1, Bk, burnin=0.3, plot=FALSE)
ss.NN <- mcmc.NN$steppingstone(10000, chain.NN, Bk, burnin=0.3, plot=FALSE)
mlnL <- c("11"=ss.11$lnr, "N1"=ss.N1$lnr, "NN"=ss.NN$lnr)
mlnL













#-------------------------------------- Concluding words
 That's probably way more than we have time for. Thanks for bearing with me. Try applying bayou to your own data, or if you don't have any, here is a script to load in a fun dataset with lots of
 possibilties. But before you finish, I have a few "words of wisdom" regarding your analyses. 

 **Caveat 1**: The reversible-jump analysis doesn't always work. Sometimes it doesn't converge in a reasonable amount of time, sometimes it doesn't find anything interesting. It's fairly
 prior sensitive to the number of shifts. Sometimes it's hard to come up with reasonable priors.
  
    + Don't put too much stock in the number of shifts recovered. Instead, consider the **number of highly supported shifts**. Consider magnitude. This means you have to understand the units
	 of the parameters. If you get about 20 shifts in your posterior and your prior was about 20 shifts, this doesn't mean you have strong support for 20 shifts. Often, no single branch will
	 have high support, and the shifts will be of negligible magnitude. However, if 5/20 DO have high support, and are of high magnitude these are likely important features to understand in your data!
  
    + In general, we found it's better to put higher priors on the number of shifts that fit models more complex than needed, than it is to put a lower prior on the number of shifts that keeps the
	 fit from becoming as complex that as it needs to be. If the prior disallows the number of shifts it needs to explain the data, often the model will collapse to a BM like model. When this happens,
	 the number of shifts and the value of theta don't really matter anymore, they don't affect anything in the model, and the model is unlikely to find the true posterior (even the *true* arrangement
	 of shifts doesn't really affect the likelihood when *alpha* is stuck in the land of Brownian Motion.) It can be useful to try many starting points, or even start with an overfitted model and let
	 *bayou* drop shifts rather than add them. 
  
    + Don't be afraid to use informative priors. Even very weakly informative priors. For example, if you're studying mammal body size, can you rule out animals larger than earth? Smaller than an atom?
	 Your prior should reflect this! You'd be surprised how often priors we use, or estimates from ML optimizations, result in parameters in these ranges. 
  
    + Think in terms of phylogenetic half-life and stationary variance. It's more intuitive. 
  
* **Caveat 2**: The reversible-jump approach is data snooping. Plus, it's hard to summarize, it's hard to conceptualize, and it's complicated. 

    + I like to take the following strategy in my analyses. I view the reversible-jump approach as an exploratory, natural history based approach for understanding the *major features* of my data.
	 I want to know which groups are unique and doing something different because we don't see in *phylogenetic time*, we see only extant species and our view of the evolutionary pattern is warped
	 by the phylogeny. rjMCMC is a *tool to see the pattern of evolution through phylogenetic time*. 
  
    + Once you identify the major features of your data, it's reasonable to ask **why?** Compare different a priori or fixed hypotheses using Bayes Factors. Find the best one as you would normally in
	 model selection. However, if that model fails to explain those *major features* of your data, you have more work to do to find a better model. 
  
    + A useful approach is to compare a model with the shifts found in your reversible-jump analysis to a model without shifts, but with an explanatory predictor added. Essentially, you ask *is the
	 data better explained by clade-level shifts (with many parameters) or a single or handful of predictor variables?* Can you *explain away* clade-level shifts with predictors? Can I *explain away*
	 several clade-specific shifts to higher regimes in whales, pinnipeds and Sirenia by including the predictor *aquatic* into my analysis?
  
* **Caveat 3**: **bayou** is not the best OU tool for every situation. **bayou** is cool because it's Bayesian. OU models are cool because they are supposed to represent biologically realistic processes
	 of adaptation. So the parameters mean something biological. So the parameters have prior expectations. Also, **bayou** is useful for exploring and understanding your data. However:

      + Right now, **bayou** doesn't fit OUwie-type models effectively. If you have questions about changing dynamics of regimes (changing *alpha* or *sig2*), such as a hypothesis that one clade is more
	 constrained than another -- use OUwie. 
  
      + There can be identifiability issues with **bayou**. You can imagine that the a similar set of shifts can result in the same distribution of data. If multiple shifts occur on a branch, for example.
	 Or if two shifts occur on neighboring branches. There are many configurations that can lead to identical or nearly identical likelihoods and mixing can be difficult. Run multiple chains. If they
	 each get stuck with alternative, but similar configurations of shifts with poor mixing, you will know that identifiability is an issue. Likelihood-based tools such as **l1ou** have implemented
	 approaches for dealing with this issue that hopefully will soon be implemented in **bayou** as well. 
  
      + The packages **slouch** and **mvSLOUCH** implement a full suite of "evolutionarily-aware" regression approaches that are likely more better than the approaches I have outlined here. They don't
	 find shifts, you have to assume them, and the can be difficult to use...but in most cases the models are more evolutionarily reasonable. These are sorely underused packages. They solve a problem
	 with PGLS that most people don't even realize is there. In **bayou**'s simple allometric models, the value of a predictor (e.g. body size) immediately effects the trait. This is a lot like PGLS
	 and reasonable if the two traits are mechanistically or developmentally linked (body mass automatically increases as you grow longer). But this isn't true for a trait like *precipitation*.
	 How much it rained this year isn't going to be a perfect predictor of body mass, though body mass may respond to precipitation. But it doesn't respond quickly. Instead, it responds to a moving
	 average of precipitation. You can't predict "optimal" body mass for a lineage from current conditions, you need the whole history of precipitation conditions to model where the optima was in the
	 past, as well as today. This means you need to model the evolutionary history of precipitation for a lineage as well as the evolutionary history of body mass. You need **slouch** to do this.
	 (Again, hopefully these models will soon be in **bayou** too)
  
      + **bayou** is not truly multivariate. Better packages for fitting multivariate OU models include **mvSLOUCH**, **mvMORPH** and **ouch**. This will also hopefully change in the future.  






























