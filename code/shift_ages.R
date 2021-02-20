library(phytools)
library(tidyverse)
library(BioGeoBEARS)

# tree simulation
str <- pbtree(n = 10, scale = 1)
plot(str, no.margin = T, show.tip.label = F)
tiplabels(col = "black", frame = "none", cex = .7, adj = c(-.5, -.5))
edgelabels(col = "red", frame = "none", cex = .7, adj = c(.5, -.5))
edgelabels(round(str$edge.length, 2), col = "darkgreen", 
           frame = "none", cex = .7, adj = c(1, 1.5))
nodelabels(col = "blue", frame = "none", cex= .7, adj = c(-.5, .5))

tr <- read.tree("data/zanne_tree_pr.tre")
branching.times(tr) %>% max()


# function to extract shift ages from branch numbers
get_shift_ages=function(tree, sh_br){
  n_hei <- nodeHeights(tree)
  shift_ages <- c()
  for(i in 1:length(sh_br)){
    e <- sh_br[i]
    shift_ages[i] <- n_hei[e,][1]+(n_hei[e,][2]-n_hei[e,][1])/2
  }
  shift_ages
}

sh_br <- c(3, 9)
get_shift_ages(str, sh_br)

# Bayou ####
get_shift_ages(shiftsum$tree, shiftsum$pars$sb) %>% round(2) %>% 
  write.csv("ou_free_fix/model_free/trait_shift_ages.csv")

# BioGeoBEARS ####
# anagenetic event time tables
ana_events_tables[[1]]$event_txt %>% head()
ana_events_tables[[1]]$event_txt %>% unique() %>% sort()

# table dimensions
tt <- list()
for(i in 1:length(ana_events_tables)){
  
  tt[[i]] <- ana_events_tables[[i]] %>% group_by(event_txt) %>% 
    summarize(mean(abs_event_time)) %>% dim()
}
tabs_list <- ana_events_tables[sapply(tt, "[", 1) == 28]

# extract average event times per dispersal type
e_ts <- tabs_list[[1]] %>% group_by(event_txt) %>% 
  summarize(mean(abs_event_time)) %>% select(event_txt) %>% data.frame()

for(i in 1:length(tabs_list)){
  
  c_temp <- tabs_list[[i]] %>% group_by(event_txt) %>% 
    summarize(mean(abs_event_time)) %>% select(`mean(abs_event_time)`) %>% 
    rename(n = `mean(abs_event_time)`)
  
  e_ts <- cbind(e_ts, c_temp)
}
names(e_ts) <- c("event_txt", paste0(rep("sm", 99), 1:99))
e_ts <- e_ts %>% mutate(mean  = round(rowMeans(across(where(is.numeric))),2) ) 
e_ts %>% write.csv("biogeob/bsm/disp_event_t.csv")

# extract mean number of events per dispersal type
f_ts <- tabs_list[[1]] %>% group_by(event_txt) %>% count() %>% 
  select(event_txt) %>% data.frame()

for(i in 1:length(tabs_list)){
  
  c_temp <- tabs_list[[i]] %>% group_by(event_txt) %>% count() %>% pull(n)
  
  f_ts <- cbind(f_ts, c_temp)
}
names(f_ts) <- c("event_txt", paste0(rep("sm", 99), 1:99))
rowMeans(f_ts)
f_ts <- f_ts %>% mutate(mean  = round(rowMeans(across(where(is.numeric))),2) ) %>% 
  select(event_txt, mean) %>% arrange(desc(mean))
f_ts %>% write.csv("biogeob/bsm/disp_event_ft.csv")

tabs_list[[1]] %>% filter(event_txt == "L->LN")

