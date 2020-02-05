
# Bayou ####
library(dplyr)
library(treeplyr)
library(bayou)
library(ggplot2)

# Model comparison
# One can fit a full model, all island together and if the posterior is over zero means that you can not regect a no effect
path <- ("custom_models/")

d <- list.files(path, all.files = T, include.dirs = T, recursive = T)
files <- d[grep("ss.", d)]
models <- unique(sapply(strsplit(files, "/"), "[[", 1))

mar_Lik <- c()
for(i in 1:length(models)){
  file <- list.files(paste0(path, models[i]), pattern = "*ss.", all.files = T)
  mar_Lik[i] <- readRDS(paste0(path, models[i],"/", file))$lnr
}

mod_comp <- data.frame(models[1:length(mar_Lik)], mar_Lik)
mod_comp$BF <-round(abs(2*(mod_comp$mar_Lik[which(max(mod_comp$mar_Lik)==mar_Lik)]-mod_comp$mar_Lik)), 2)
mod_comp[which(max(mod_comp$BF)==mod_comp$BF),] # Best model
best_model <- mod_comp[which(max(mod_comp$BF)==mod_comp$BF),][1,1]
#

chain.N1111 <- readRDS(paste0(path, best_model, "/chain.", best_model, ".rds"))
chain.N1111 <- set.burnin(chain.N1111, 0.3)
sum_N1111 <- summary(chain.N1111)
sum_N1111 <- readRDS(paste0(path, best_model, "/sum_", best_model, ".rds"))
sum_N1111$statistics
head(sum_N1111$branch.posteriors)

plot(chain.N1111, auto.layout = FALSE)

# regression values
par_names <- rownames(sum_N1111$statistics)[grep("beta_", rownames(sum_N1111$statistics))]
x_names <- sapply(strsplit(par_names, "_"), "[[", 2)
tab_sum <- data.frame(x = x_names,
           beta = sum_N1111$statistics[par_names,"Mean"],
           sd = sum_N1111$statistics[par_names,"SD"],
           hpdL = sum_N1111$statistics[par_names,"HPD95Lower"],
           hpdU = sum_N1111$statistics[par_names,"HPD95Upper"])

png("res_model.png", width = 900, height = 900, bg = "transparent", res = 250)
ggplot(tab_sum, aes(x = x, y = beta)) + geom_point() +
  geom_errorbar(aes(ymin = hpdL, ymax = hpdU), width = .2) + theme_classic() + theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(x = "", y = "Island Effect (β)")
dev.off()

plotBranchHeatMap(d_fruit_lg_ln$phy, chain.N1111, "alpha", pal = cm.colors, cex = .1)
plotBranchHeatMap(d_fruit_lg_ln$phy, chain.N1111, "beta_newguinea", pal = cm.colors, cex = .1)

# l1ou ####
library(dplyr)
library(treeplyr)
library(genlasso)
library(l1ou)
library(geiger)
library(phytools)

e.Model<-readRDS("nl1ou/output/eModel.rds")
eModel_boot_sup<-readRDS("nl1ou/output/eModel_boot_sup.rds")
eModel_c<-readRDS("nl1ou/output/eModel_c.rds")
eMode_cd<-readRDS("nl1ou/output/eModel_cd.rds")
