library(tidyverse)

med_tetha <- sapply(shiftsum$cladesummaries, unlist)[1,] %>% sapply("[") #106
sb_shifts <- c(0, shiftsum$pars$sb) # 106
shiftsum$descendents # 106

cons_tab <- data.frame()
for(i in 1:length(med_tetha)){
  temp <- data.frame(spp = shiftsum$descendents[[i]], 
             median_theta = rep(med_tetha[i], length(shiftsum$descendents[[i]])), 
             shift_sb = rep(sb_shifts[i], length(shiftsum$descendents[[i]])))
  cons_tab <- rbind(cons_tab,temp)
}

cons_tab %>% tail()

tip_path <- get_path(1:Ntip(shiftsum$tree), shiftsum$tree$edge, Ntip(shiftsum$tree)+1)

n_shifts <- c()
for(i in 1:length(tip_path)){
  n_shifts[i] <- table(sb_shifts %in% tip_path[[i]])["TRUE"]
}
n_shifts[is.na(n_shifts)] <- 0

data.frame(spp = shiftsum$tree$tip.label, n_shifts)

res_cons_tab <- left_join(cons_tab, data.frame(spp = shiftsum$tree$tip.label, n_shifts))
write.csv(res_cons_tab, file = "sum_shifts_indom_.csv")
res_cons_tab %>% tail()

###
par(mfrow=c(1,2))
plot((shiftsum$tree), main = "no- lad", show.tip.label = F, no.margin = T)
plot(ladderize(shiftsum$tree), main = "lad", show.tip.label = F, no.margin = T)

tt <- phytools::pbtree(n = 10, scale = 1)
plot((tt))
tiplabels(col = "blue", frame = "none", offset = .1)
edgelabels(col = "red", frame = "none", cex = .7)
get_path(1:Ntip(tt), tt$edge, Ntip(tt)+1)
