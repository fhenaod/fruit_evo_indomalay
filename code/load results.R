
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

mod_comp <- data.frame(models = models[1:length(mar_Lik)], mar_Lik)
mod_comp$BF <-round(abs(2*(mod_comp$mar_Lik[which(max(mod_comp$mar_Lik, na.rm = T)==mar_Lik)]-mod_comp$mar_Lik)), 2)
mod_comp[which(max(mod_comp$BF, na.rm = T)==mod_comp$BF),] # Best model
mod_comp %>% arrange((BF))
best_model <- mod_comp[which(max(mod_comp$BF, na.rm = T)==mod_comp$BF),][1,1]
#

best_model <- "N1111"
chain.best <- readRDS(paste0(path, best_model, "/chain.", best_model, ".rds"))
chain.best <- set.burnin(chain.best, 0.3)
#sum_best <- summary(chain.best)
sum_best <- readRDS(paste0(path, best_model, "/sum_", best_model, ".rds"))
sum_best$statistics
head(sum_best$branch.posteriors)

plot(chain.best, auto.layout = FALSE)

# regression values
par_names <- rownames(sum_best$statistics)[grep("beta_", rownames(sum_best$statistics))]
x_names <- sapply(strsplit(par_names, "_"), "[[", 2)
tab_sum <- data.frame(x = x_names,
           beta = sum_best$statistics[par_names,"Mean"],
           sd = sum_best$statistics[par_names,"SD"],
           hpdL = sum_best$statistics[par_names,"HPD95Lower"],
           hpdU = sum_best$statistics[par_names,"HPD95Upper"])

png("res_model.png", width = 900, height = 900, bg = "transparent", res = 250)
ggplot(tab_sum, aes(x = x, y = beta)) + geom_point() +
  geom_errorbar(aes(ymin = hpdL, ymax = hpdU), width = .2) + theme_classic() + theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(subtitle = paste("Model:",best_model), x = "", y = "Island Effect (β)")
dev.off()

plotBranchHeatMap(d_fruit_lg_ln$phy, chain.best, "alpha", pal = cm.colors, cex = .1)
plotBranchHeatMap(d_fruit_lg_ln$phy, chain.best, "beta_newguinea", pal = cm.colors, cex = .1)
plotBranchHeatMap(d_fruit_lg_ln$phy, chain.best, "theta", pal = cm.colors, cex = .1)

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
