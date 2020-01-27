
# Bayou ####
library(dplyr)
library(treeplyr)
library(bayou)
library(ggplot2)

# Model comparison
# One can fit a full model, all island together and if the posterior is over zero means that you can not regect a no effect
models <- c("N1000", "N1100", "N1110", "N1111")
path <- ("fixed_ou/")
mar_Lik <- c()
for(i in 1:length(models)){
  file <- list.files(paste0(path, models[i]), pattern = "*ss.N", all.files = T)
  mar_Lik[i] <- readRDS(paste0(path, models[i],"/", file))$lnr
}

mod_comp <- data.frame(models, mar_Lik)
mod_comp$BF <-round(abs(2*(-1899.054-mod_comp$mar_Lik)), 2)
best_mod <- "N1111"

#
chain.N1111 <- readRDS("fixed_ou/N1111/chain.N1111.rds")
chain.N1111 <- set.burnin(chain.N1111, 0.3)
sum_N1111 <- summary(chain.N1111)
sum_N1111 <- readRDS("fixed_ou/N1111/sum_N1111.rds")
sum_N1111$statistics

plot(chain.N1111, auto.layout = FALSE)

# regression values
tab_sum <- data.frame(x = c("sunda", "sulawesi", "malaku", "newguinea"),
           beta = sum_N1111$statistics[c("beta_sunda", "beta_sulawesi", "beta_maluku", "beta_newguinea"),"Mean"],
           sd = sum_N1111$statistics[c("beta_sunda", "beta_sulawesi", "beta_maluku", "beta_newguinea"),"SD"],
           hpdL = sum_N1111$statistics[c("beta_sunda", "beta_sulawesi", "beta_maluku", "beta_newguinea"),"HPD95Lower"],
           hpdU = sum_N1111$statistics[c("beta_sunda", "beta_sulawesi", "beta_maluku", "beta_newguinea"),"HPD95Upper"])

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
