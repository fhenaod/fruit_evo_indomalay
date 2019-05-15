# Modifed from Lizzie Wolkovich's github code for Davies etal (2018; 
# 'Phylogenetically weighted regression'); https://github.com/lizzieinvancouver/pwr



suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(nlme))
suppressPackageStartupMessages(require(phylobase))
suppressPackageStartupMessages(require(grid))
suppressPackageStartupMessages(require(subplex))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(phytools))
suppressPackageStartupMessages(require(dplyr))



# Desktop
#setwd('D:/Box Sync/Projects/Projects (active)/Indo archipelago survey/Flora Malesiana analysis/MS Tree fruit evolution/Analysis')
# Laptop
setwd('C:/Users/JB/Box Sync/Projects/Projects (active)/Indo archipelago survey/Flora Malesiana analysis/MS Tree fruit evolution/Analysis')



#--------- DATA AND FUNCTIONS -----------------------------
rm(list = ls())
source("0pwr-functions.R")
source("0pwr-plots.R")
spptree1 <- ape::read.tree("spptreeZanneShort.tree") 
spptraits1 <- read.csv("spptraits1.csv")
alldata <- read.csv("FloraMalesiana190226.csv")
superfamilies <- data.frame(read.csv("superfamilies.csv"))

# clean and link data
spptraits <- cbind.data.frame(family=spptraits1$family, genus=spptraits1$genus, Gen_spp=spptraits1$Gen_spp, 
	sunda=spptraits1$sunda, sulawesi=spptraits1$sulawesi, maluku=spptraits1$maluku, newguinea=spptraits1$newguinea,
	Y=spptraits1$fruit_lg)
spptraits <- spptraits[complete.cases(spptraits),]
colnames(alldata)[colnames(alldata)=="Family"] <- "family"
alldata$family <- as.character(alldata$family)
spptraits$family <- as.character(spptraits$family)
superfamilies$family <- as.character(superfamilies$family)
alldata$family <- as.character(alldata$family)
spptraits$family <- as.character(spptraits$family)
spptraits <- dplyr::left_join(spptraits, superfamilies, by="family")
alldata <- dplyr::left_join(alldata, superfamilies, by="family")

# Remove species in trait dataset that are not in the phylogeny
missing.i <- which(!spptraits$Gen_spp %in% spptree1$tip.label)
spptraits <- spptraits[-missing.i, ]

# Remove species in phylogeny that are not in trait data
spptree <- spptree1
in.spec <- spptree$tip.label %in% spptraits$Gen_spp
spptree <- drop.tip(spptree, spptree$tip.label[!in.spec])




#---------------- FRUIT LENGTH ANALYSIS ------------------------------------------------------------
# merge data to generate phylo4d object
datp4d <- phylo4d(spptree, spptraits, label.type="column", label.column="Gen_spp")

#--------- PGLS 
pgls.b <- gls(Y ~ sunda + sulawesi + maluku + newguinea, data=tipData(datp4d), correlation=corBrownian(phy=spptree)) # brownian model
pgls.m <- gls(Y ~ sunda + sulawesi + maluku + newguinea, data=tipData(datp4d), correlation=corMartins(value=1, phy=spptree)) # O-U model
pgls.b
pgls.m
coef(pgls.m)
confint(pgls.m)

#--------- PWR (uses functions defined by source: "0pwr-functions.R") 
# Can't do model selection on these PWR models. Instead, calculate the PWR version of the "best" PGLS model, above
pwr.b <- pwr(Y ~ sunda + sulawesi + maluku + newguinea, datp4d, wfun="brownian")
#pwr.g <- pwr(Y ~ sunda + sulawesi + maluku + newguinea, datp4d, wfun="gaussian")
# Get optimal bandwidth (measured in branch-length units) for O-U model
#bw <- get.opt.bw(Y ~ sunda + sulawesi + maluku + newguinea, datp4d, wfun="martins", method="subplex")
bw <- 5.313353
pwr.m <- pwr(Y ~ sunda + sulawesi + maluku + newguinea, datp4d, bw=bw, wfun="martins")

#--------- Models with 1 variable at a time
pgls1 <- gls(Y ~ sunda, data=tipData(datp4d), correlation=corMartins(value=1, phy=spptree))
pgls2 <- gls(Y ~ sulawesi, data=tipData(datp4d), correlation=corMartins(value=1, phy=spptree))
pgls3 <- gls(Y ~ maluku, data=tipData(datp4d), correlation=corMartins(value=1, phy=spptree))
pgls4 <- gls(Y ~ newguinea, data=tipData(datp4d), correlation=corMartins(value=1, phy=spptree))
pwr1 <- pwr(Y ~ sunda, datp4d, bw=bw, wfun="martins")
pwr2 <- pwr(Y ~ sulawesi, datp4d, bw=bw, wfun="martins")
pwr3 <- pwr(Y ~ maluku, datp4d, bw=bw, wfun="martins")
pwr4 <- pwr(Y ~ newguinea, datp4d, bw=bw, wfun="martins")
 
# Coefficient estimates
# est=regression coefficient for each spp, lb=lower bound, ub=upper bound
pwrests1 <- getEst(pwr1)
pwrests2 <- getEst(pwr2)
pwrests3 <- getEst(pwr3)
pwrests4 <- getEst(pwr4)
summary(pwrests1) 
summary(pwrests2) 
summary(pwrests3) 
summary(pwrests4)

#--------- Correlations
cortest <- cbind(sunda=pwrests1$est, sula=pwrests2$est, maluku=pwrests3$est, ng=pwrests4$est)
cor(cortest)

#--------- TREE FIGURES (uses functions defined by source: "0pwr-plots.R") 
modpgls <- pgls1
modpwr <- pwr1
pdf("treefig_test.pdf", width=8, height=50) # width & height in inches
tp(addData(datp4d[,0], data.frame(
        gest = coef(modpgls)[2],
        glb = confint(modpgls)[2,1],
        gub = confint(modpgls)[2,2],
        getEst(modpwr)    )),
    show.estimates=TRUE, aex=1.2)
dev.off()






#--------- Are younger or more speciose families more evol labile?
datfam <- data.frame(family=unique(spptraits$family), numsppFM=NA, numtips=NA, nodenum=NA, age=NA, coef.mean=NA, coef.se=NA)
datfam$family <- as.character(datfam$family)

#modpwr <- pwr1
#modpwr <- pwr2
#modpwr <- pwr3
#modpwr <- pwr4

sppcoefs <- cbind.data.frame(Gen_spp=spptree$tip.label, coef=getEst(modpwr)$est)
sppcoefs$Gen_spp <- as.character(sppcoefs$Gen_spp)
tmp <- suppressWarnings(dplyr::left_join(sppcoefs, spptraits, by="Gen_spp"))
sppcoefs <- subset(tmp, select=-c(genus, sunda, sulawesi, maluku, newguinea, Y))

for(i in 1:nrow(datfam)){
   tryCatch({
	#i=4 
	datfam$numsppFM[i] <- nrow(subset(alldata, family==datfam$family[i]))# num spp in the family (in Flora Malesiana)
	famspp <- subset(spptraits, family==datfam$family[i]) # subset of traits database for species in family i
	datfam$nodenum[i] <- ape::getMRCA(spptree, tip=c(as.character(famspp$Gen_spp))) # node number of most recent ancestor of family i
	famtree <- extract.clade(spptree, datfam$nodenum[i], root.edge=0) # phylogeny of family i
	datfam$numtips[i] <- length(famtree$tip.label) # number of tips in the phylogeny of family i
	datfam$age[i] <- max(branching.times(famtree)) # age of family i (in branch-length units)	
	sppcoefs.i <- subset(sppcoefs, family==datfam$family[i]) # regression coefficients for family i
	datfam$coef.mean[i] <- mean(sppcoefs.i$coef) # mean regression coefficient for the family	
	datfam$coef.se[i] <- sd(sppcoefs.i$coef)/sqrt(nrow(sppcoefs.i)) # SE of regression coef for the family
   }, error=function(e){}) # end of the tryCatch function
}
datfam <- datfam[complete.cases(datfam),]

# quality check
datfam$bad <- ifelse(datfam$numtips > datfam$numsppFM, 1, 0)
BadFams <- subset(datfam, bad ==1)[,1]
datfam <- subset(datfam, bad != 1) # identifies Flora Malesiana families that are polyphyletic

# some plots
plot(datfam$numsppFM, datfam$coef.mean)
plot(log10(datfam$numsppFM), datfam$coef.mean)
plot(datfam$age, datfam$coef.mean)

# some stats
summary(lm(coef.mean ~ log10(numsppFM), data=datfam))





#--------- Are younger or more speciose superfamilies more evol labile?
datsup <- data.frame(superfamily=unique(superfamilies$superfamily), numsppFM=NA, numtips=NA, nodenum=NA, age=NA, coef.mean=NA, coef.se=NA)
datsup$superfamily <- as.character(datsup$superfamily)

# remove bad families
for(j in 1:length(BadFams)){spptraits <- subset(spptraits, family != BadFams[j])	}

#modpwr <- pwr1 # sunda
#modpwr <- pwr2 # sulawesi
#modpwr <- pwr3 # maluku
#modpwr <- pwr4 # new guinea

sppcoefs <- cbind.data.frame(Gen_spp=spptree$tip.label, coef=getEst(modpwr)$est)
sppcoefs$Gen_spp <- as.character(sppcoefs$Gen_spp)
tmp <- suppressWarnings(dplyr::left_join(sppcoefs, spptraits, by="Gen_spp"))
sppcoefs <- subset(tmp, select = -c(family, genus, sunda, sulawesi, maluku, newguinea, Y))

for(i in 1:nrow(datsup)){
   tryCatch({
	#i=1
	datsup$numsppFM[i] <- nrow(subset(alldata, superfamily==datsup$superfamily[i]))# num spp in the superfamily (in Flora Malesiana)
	superfamspp <- subset(spptraits, superfamily==datsup$superfamily[i]) # subset of traits database for species in superfamily i
	datsup$nodenum[i] <- ape::getMRCA(spptree, tip=c(as.character(superfamspp$Gen_spp))) # node number of most recent ancestor of family i
	superfamtree <- extract.clade(spptree, datsup$nodenum[i], root.edge=0) # phylogeny of superfamily i
	datsup$numtips[i] <- length(superfamtree$tip.label) # number of tips in the phylogeny of superfamily i
	datsup$age[i] <- max(branching.times(superfamtree)) # age of superfamily i (in branch-length units)	
	sppcoefs.i <- subset(sppcoefs, superfamily==datsup$superfamily[i]) # regression coefficients for superfamily i
	datsup$coef.mean[i] <- mean(sppcoefs.i$coef) # mean regression coefficient for the superfamily	
	datsup$coef.se[i] <- sd(sppcoefs.i$coef)/sqrt(nrow(sppcoefs.i)) # SE of regression coef for the superfamily
   }, error=function(e){}) # end of the tryCatch function
}
datsup <- datsup[complete.cases(datsup),]

# quality check
datsup$bad <- ifelse(datsup$numtips > datsup$numsppFM, 1, 0)
datsup <- subset(datsup, bad != 1) # identifies Flora Malesiana superfamilies that are polyphyletic

# Redo the clade ages (MYA; from Magallon et al. 2015; New Phytologist; Table S2; 'Crown' ages)
a1 <- 122.58 # Superasteridae # 119.55 125.55
a2 <- 122.41 # Superrosidae # 119.49 125.38
a3 <- 114.83 # Ranunculales # 112.14 123.18
a4 <- 122.56 # Chloranthales # 120.84 126.74
a5 <- NA
a6 <- 132.39 # Magnoliidae # 130.24 134.14
a7 <- 117.38 # Proteales 109.56 125.73
a8 <- 133.23 # Monocots # 131.61 134.69
datsup$age_manual <- c(a1, a2, a3, a4, a5, a6, a7, a8)

# some stats
summary(lm(coef.mean ~ numsppFM, data=datsup))
summary(lm(coef.mean ~ log10(numsppFM), data=datsup))
summary(lm(coef.mean ~ age_manual, data=datsup))
summary(lm(coef.mean ~ log10(age_manual), data=datsup))
summary(lm(coef.mean ~ age, data=datsup))
summary(lm(coef.mean ~ log10(age), data=datsup))

# some plots
plot(datsup$numsppFM, datsup$coef.mean)
plot(log10(datsup$numsppFM), datsup$coef.mean)
plot(datsup$age_manual, datsup$coef.mean)
plot(log10(datsup$age_manual), datsup$coef.mean)

# ANOVA of major lineages
sppcoefs$superfamily <- as.factor(sppcoefs$superfamily)
datshort <- subset(sppcoefs, superfamily != "")
modaov <- aov(coef ~ superfamily, data=datshort)
summary(modaov)
TukeyHSD(modaov)


# Figures & Tables
# TukeyHSD(modaov) exported to GraphsTables.xls>"TukeyHSD"
# datsup exported to GraphsTables.xls>"Bar Graphs"












#--------- Phylogeny plots
plot.phylo(spptree, show.tip.label=F)
nodelabels(, col="black", bg="gray")










#--- interpreting the bandwidth parameter (from Davies etal, p.9)
# "We suggest that this bandwidth might be broadly interpreted as indicative of the evolutionary rate of change
# in model coefficients and can be likened to the phylogenetic half-life concept of Hansen (1997). A narrow bandwidth may reflect faster
# rates of change, such that modelled evolutionary relationships are divergent even among closely related species. Conversely, a broad
# bandwidth would indicate slower rates of change"








#--------- SCATTERPLOT OF DATA (uses functions defined by source: "0pwr-plots.R") -----------------------
# See Fig 3 and Section 2.2.1 of Davies etal for an explanation of this plot
# The open & closed circles refer to the two clades after the most ancient bifurcating split

# points labelled by species name
pdf("scatter.pdf", width=5, height=5)
par(xpd=TRUE)
plot(tipData(datp4d)$seed, tipData(datp4d)$FFD, pch=rep(c(1,16), each=4), 
      xlab="Predictor", ylab="Response", bty="l")
	text(tipData(datp4d)$seed, tipData(datp4d)$FFD, labels=rownames(tipData(datp4d)), pos=3, cex=0.8)
dev.off()

# points unlabelled
pdf("scatter.pdf", width=5, height=5)
par(xpd=TRUE)
plot(tipData(datp4d)$seed, tipData(datp4d)$FFD, pch=rep(c(1,16), each=4), 
      xlab="Predictor", ylab="Response", bty="l")
	#text(tipData(datp4d)$seed, tipData(datp4d)$FFD, labels=rownames(tipData(datp4d)), pos=3, cex=0.8)
dev.off()































