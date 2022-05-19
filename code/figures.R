library(tidyverse)
library(bayou)
library(BioGeoBEARS)
library(castor)

chain_fixed <- readRDS("run2/ou_free_fix/model_free/mcmc_fixed.pp75.rds")

# Fig 2 ####
data_spp <- read.csv(file = "data/species_clades.csv", head = T)
data_spp %>% head()
seed_tax <- read.csv(file = "data/seed_tax.csv", head = T)

left_join(
data_spp %>% 
  mutate(genus = sapply(strsplit(data_spp$species, "_", fixed = T), "[", 1)),
seed_tax, by = c("genus" = "genus")
) %>% filter(clade == "Fabids") %>% filter(order == "Huerteales")
  pull(order) %>% unique()
  filter(order == "Malpighiales" | order == "Fabales") 

data_spp <- 
data_spp %>% 
  mutate(clade = ifelse(species == "Irvingia_malayana" |
                        species == "Suriana_maritima", "Fabids", clade),
         clade = ifelse(species == "Perrottetia_alpestris" |
                        species == "Perrottetia_moluccana", "Malvids", clade))

par(mfrow = c(1,2)) 
plotSimmap.mcmc(chain_fixed, edge.type = "theta",
                pp.cutoff = .3, cex = .001, no.margin = F, circles = F,
                pal = terrain.colors,
                legend_settings = 
                  list(plot = T, x = 0.1, y = 500, 
                       cex.lab = .5, adjx = .02))
mtext(LETTERS[1], side = 3, adj = .05, line = -1.3, cex = 1.5, font = 2)
axisPhylo()

tree2paint <- d_fruit_lg_ln$phy
clad_nod <- data_spp %>% pull(clade) %>% unique()
clad_nod <- clad_nod[c(1:3,6)]
for(i in 1:length(clad_nod)){
  clad2name <- data_spp %>% filter(clade == clad_nod[i]) %>% 
    select(species) %>% c()
  spp_mrca <- castor::get_mrca_of_set(tree2paint, clad2name$species)
  print(paste0("node mrca ", clad_nod[i],": " , spp_mrca))
  tree2paint <- paintSubTree(tree2paint, spp_mrca, 
               state = clad_nod[i], anc.state = "0", stem = T)
}

cols <- setNames(c("black", viridis::viridis(length(clad_nod))), c("0", clad_nod))
png(file ="fig2_a.png", 
    width = 25, height = 25, res = 300, units = "cm")
plot(tree2paint, fsize = 0.001, 
     colors = cols, direction = "rightwards")

add.simmap.legend(colors = setNames(cols[2:length(cols)], clad_nod), 
                  fsize = 0.35, prompt = FALSE, x = 5, y = 300)
mtext(LETTERS[2], side = 3, adj = .05, line = -1.3, cex = 1.5, font = 2)
dev.off()

nodelabels("F", node =  3079, frame = "none")
nodelabels("M", node =  3894, frame = "none")

# if tree in fan type
par(mfrow = c(1, 1))
plotSimmap.mcmc(chain_free, edge.type = "theta", type = "fan",
                pp.cutoff = .3, cex = .01, no.margin = T, circles = F,
                pal = terrain.colors, 
                legend_settings = 
                  list(plot = T, x = -1.1, y = .4, 
                       cex.lab = .5, adjx = .02))

sup_clad_nod <- unique(data_spp$superclade)
# clade names in nodes 
for(i in 1:length(sup_clad_nod)){
  clad2name <- data_spp %>% filter(superclade == sup_clad_nod[i]) %>% 
    select(species) %>% c()
  spp_mrca <- castor::get_mrca_of_set(tree2paint, clad2name$species)
  print(paste0("node mrca ", sup_clad_nod[i],": " , spp_mrca))
  nodelabels(sup_clad_nod[i], spp_mrca,
             #adj = c(1.1, 0),
             frame = "n", cex = 1, font = 3)
}  

# clade names in lines arround tips
cols <- setNames(c(viridis::viridis(length(sup_clad_nod))), c(sup_clad_nod))
cols <- setNames(c(RColorBrewer::brewer.pal(n = length(sup_clad_nod), name = "Set1")), c(sup_clad_nod))

for(i in 1:length(sup_clad_nod)){
  clad2name <- data_spp %>% filter(superclade == sup_clad_nod[i]) %>% 
    select(species) %>% c()
  
  arc.cladelabels(d_fruit_lg_ln$phy, sup_clad_nod[i], 
                  get_mrca_of_set(d_fruit_lg_ln$phy, clad2name$species),
             orientation = "horizontal", mark.node = T,
             ln.offset = 1.02, lab.offset = 1.04, cex = .9, lwd = 4,
             wing.length = 0.001, col = cols[i])
} 

# using ggtree
library(ggtree)

data_spp %>% head()
sup_clad <- data_spp$superclade %>% unique()

nodes2paint <- c()
for(i in 1:length(sup_clad)){
  nodes2paint[i] <- 
    data_spp %>% 
    filter(superclade == sup_clad[i]) %>% pull(species) %>% 
    castor::get_mrca_of_set(tree2paint, .)
}
nodes_dt <- data.frame(sup_clad, nodes2paint)

tree2paint %>% groupClade(., .node = nodes2paint[c(1,2,4)]) %>% 
  ggtree(aes(color = group), ladderize = F) +
  #theme_tree2() +

ggtree(tree2paint, ladderize = F) +
  geom_balance(node = nodes2paint[1], fill = 'steelblue', color = 'white', alpha = 0.2, extend = 1) +
  geom_balance(node = nodes2paint[2], fill = 'darkgreen', color = 'white', alpha = 0.2, extend = 1) +
  geom_balance(node = nodes2paint[4], fill = 'coral', color = 'white', alpha = 0.2, extend = 1) +
  geom_cladelab(node = nodes2paint[1], label = paste0(sup_clad[1]), align = TRUE, angle = 270, fontsize = 3.5, offset = .01, offset.text = .25, textcolor = 'steelblue', barcolor = 'steelblue') +
  geom_cladelab(node = nodes2paint[2], label = paste0(sup_clad[2]), align = TRUE, angle = 270, fontsize = 3.5, offset = .01, offset.text = .25, textcolor = 'darkgreen', barcolor = 'darkgreen') +
  geom_cladelab(node = nodes2paint[4], label = paste0(sup_clad[4]), align = TRUE, angle = 270, fontsize = 3.5, offset = .01, offset.text = .25, textcolor = 'coral', barcolor = 'coral') 
  # + theme_tree2()

# Fig 3 ####
tr <- d_fruit_lg_ln$phy
par(mfrow = c(1,1))
analysis_titletxt <- ""
results_object <- resBAYAREALIKE

png(file ="fig3_a.png", 
    width = 25, height = 25, res = 300, units = "cm")
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
                         addl_params = list("j"), plotwhat = "pie",
                         tipcex = 0, 
                         label.offset = 0.45, show.tip.label = F, 
                         statecex = 0.4, splitcex = 0.6, titlecex = 0.001,
                         plotsplits = F, cornercoords_loc = scriptdir, 
                         include_null_range = T, tr = tr, tipranges = tipranges, 
                         plotlegend = F, tipboxes_TF = F)
dev.off()

png(file ="fig3_b.png", 
    width = 25, height = 25, res = 300, units = "cm")
plotSimmap.mcmc(chain_fixed, edge.type = "regimes", 
                pp.cutoff = 0.9, cex = .001, no.margin = T, 
                direction = "leftwards", circles = F)
dev.off()

# transistion times histo-densi graph ####
island_times <- read.csv("run2/biogeob/bsm_bayes/disp_event_final_tab.csv", 
                         header = T)

fruit_times <- read.csv("run2/ou_free_fix/model_free/trait_shift_ages.csv", 
                        header = T)

rbind(
island_times %>% select(-sd_time_events, -mean_n_events) %>% #head()
  gather(key, value, -event_txt) %>% 
  select(value, key),
fruit_times %>% mutate(value = x, key = "fruit_change") %>% 
  select(value, key)
) %>% ggplot(aes(x = value, color = key, fill = key)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5) +
  geom_density(alpha = .6) +
  #geom_histogram(position = "dodge") + 
  #geom_vline()
  theme_classic() +
  labs(x = "Millons of years", y = "Density", )


# Fig 5. #####
# island
f5a <- 
read.csv("run2/custom_models/fixed_island_fix/tab_sum_island_fixed.csv", 
         header = T) %>% 
  ggplot(aes(x = factor(x, level = c('sunda', 'sulawesi', 'maluku','newguinea')), y = beta)) + 
  geom_point() + geom_errorbar(aes(ymin = hpdL, ymax = hpdU), width = .2) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none",
        axis.text.x = 
          element_text(angle = 0, vjust = 0.9, hjust = .5, colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(x = "", y = "β coefficient") + # Island effect (β)
  scale_x_discrete(labels = c('Sunda', 'Sulawesi', 'Maluku', 'New  \nguinea'))
ggsave("fig5_a.png", width = 15, height = 12, units = "cm")

# habit
f5b <- 
read.csv("run2/custom_models/fixed_habit_fix/tab_sum_habit_fixed.csv", 
         header = T) %>% 
  ggplot(aes(x = factor(x, level = c('climber_woody', 'climber_herbaceous', 'woody', 'herbaceous')), y = beta)) + 
  geom_point() + geom_errorbar(aes(ymin = hpdL, ymax = hpdU), width = .2) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none",
        axis.text.x = 
          element_text(angle = 0, vjust = 0.9, hjust = 0.5, colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(x = "", y = "") + # Habit effect (β)
  scale_x_discrete(labels = c('Climber\nwoody', 'Climber\nherbaceous',  'Woody', 'Herbaceous'))
ggsave("fig5_b.png", width = 15, height = 12, units = "cm")

# type
f5c <- 
read.csv("run2/custom_models/fixed_type_fix/tab_sum_type_fixed.csv", 
         header = T) %>% 
  ggplot(aes(x = factor(x, level = c('fleshy', 'dry')), y = beta)) + 
  geom_point() + geom_errorbar(aes(ymin = hpdL, ymax = hpdU), width = .5) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none",
        axis.text.x = 
          element_text(angle = 0, vjust = 0.9, hjust = .5, colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(x = "", y = "") + #Fruit type effect (β)
  scale_x_discrete(labels = c('Fleshy\n', 'Dry'))
ggsave("fig5_c.png", width = 15, height = 12, units = "cm")

ggpubr::ggarrange(f5a, f5b, f5c, nrow = 1, labels = c("A", "B", "C"))
ggsave("fig5_panel.png", width = 38, height = 12, units = "cm")
