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
path <- c("run2/ou_free_fix/model_free/")
shiftsum <- readRDS(path0(path, "shift_sum_fixed.pp75.rds"))
get_shift_ages(shiftsum$tree, shiftsum$pars$sb) %>% round(2) %>% 
  write.csv(paste0(path,"trait_shift_ages.csv"))

# BioGeoBEARS ####
# anagenetic event time tables
ana_events_tables[[1]]$event_txt %>% head()
ana_events_tables[[1]]$event_txt %>% unique() %>% sort()

# tables dimensions
tt <- list()
for(i in 1:length(ana_events_tables)){

    tt[[i]] <- ana_events_tables[[i]] %>% group_by(event_txt) %>% 
    summarize(mean(abs_event_time)) %>% dim()
}
max_n_evts <- which(sapply(tt, "[", 1)==max(sapply(tt, "[", 1)))

# extract average event times per dispersal type
for(i in 1:length(ana_events_tables)){
  if(i == 1){
    c_temp <- left_join(
        ana_events_tables[[max_n_evts]] %>% group_by(event_txt) %>% 
          summarize(mean(abs_event_time)) %>% select(event_txt) %>% data.frame(),
        ana_events_tables[[i]] %>% group_by(event_txt) %>% 
          summarize(mean(abs_event_time)) %>% rename(n = `mean(abs_event_time)`), 
        by = "event_txt")
  } else {
    c_temp <- left_join(
      c_temp,
      ana_events_tables[[i]] %>% group_by(event_txt) %>% 
        summarize(mean(abs_event_time)) %>% rename(n = `mean(abs_event_time)`), 
      by = "event_txt")
  } 
}
names(c_temp) <- c("event_txt", paste0(rep("sm", length(ana_events_tables)), 1:length(ana_events_tables)))
e_ts <- c_temp
e_ts <- e_ts %>% 
  mutate(mean_time_events = round(rowMeans(na.rm = T, across(where(is.numeric))),2)) %>% 
  group_by(event_txt) %>% mutate(sd_time_events = round(sd(na.rm = T, across(where(is.numeric))),2) )

e_ts %>% select(event_txt,mean_time_events, sd_time_events) %>% 
  write.csv("disp_event_times.csv")

# extract mean number of events per dispersal type
for(i in 1:length(ana_events_tables)){
  if(i == 1){
    n_temp <- 
      left_join(
        ana_events_tables[[max_n_evts]] %>% group_by(event_txt) %>% 
          count() %>% select(event_txt) %>% data.frame(),
        ana_events_tables[[i]] %>% group_by(event_txt) %>% count(), 
        by = "event_txt"
      )
  } else {
    n_temp <- left_join(
      n_temp,
      ana_events_tables[[i]] %>% group_by(event_txt) %>% count(), 
      by = "event_txt"
    )
  } 
}
names(n_temp) <- c("event_txt", paste0(rep("sm", length(ana_events_tables)), 1:length(ana_events_tables)))
f_ts <- n_temp
f_ts <- f_ts %>% mutate(mean_n_events  = round(rowMeans(na.rm = T, across(where(is.numeric))),2) )

path_out <- c("run2/biogeob/bsm_bayes/")
f_ts %>% select(event_txt, mean_n_events) %>% 
  write.csv(paste0(path_out, "disp_event_count.csv"))

left_join(e_ts %>% select(event_txt,mean_time_events, sd_time_events), 
          f_ts %>% select(event_txt, mean_n_events), by = "event_txt") %>% 
  write.csv(paste0(path_out, "disp_event_final_tab.csv"), row.names = F)
