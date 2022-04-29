library(tidyverse)
library(V.PhyloMaker)

spptraits1 <- read.csv("data/data210630.csv")
#spptraits1$species[grep("-_", spptraits1$species)]
df2tree <- spptraits1 %>% 
  mutate(g_dat = sunda + sulawesi + maluku + newguinea) %>% 
  filter(g_dat > 0) %>% 
  select(species = Gen_spp, genus = genus, family = family)

tree_n1_s3 <- phylo.maker(sp.list = df2tree, tree = GBOTB.extended, 
                       nodes = nodes.info.1, scenarios = "S3")

tree_n2_s3 <- phylo.maker(sp.list = df2tree, tree = GBOTB.extended, 
                          nodes = nodes.info.2, scenarios = "S3")

plot(tree_n2_s3$scenario.3, #type = "fan", 
     show.tip.label = F, no.margin = T)
nodelabels(col = "red", frame = "none", cex = .5)

spptree1 <- tree_n2_s3$scenario.3 %>% multi2di()
spptree1$edge.length[spptree1$edge.length<=0] <- 0.03
spptree1$edge.length <- spptree1$edge.length + 0.03
#spptree1 <- phytools::force.ultrametric(spptree1, method = "extend")
spptree1 <- castor:::extend_tree_to_height(spptree1)$tree

table(spptree1$edge.length<=0)
min(spptree1$edge.length)

spptree1 %>% is.ultrametric(tolerance = 0) ; spptree1 %>% is.binary() ; spptree1 %>% is.rooted()
write.tree(spptree1, file = "data/indomalay_s3_v2.tre")

## #####
sp_tre_raw <- read.tree("data/spptreeSmithBrownShort.tre") 
min(sp_tre_raw$edge.length)
sp_tre_raw$edge.length[sp_tre_raw$edge.length<=0] <- 0.01
sp_tre_raw <- castor:::extend_tree_to_height(sp_tre_raw)$tree

table(sp_tre_raw$edge.length<=0)
min(sp_tre_raw$edge.length)
write.tree(sp_tre_raw, file = "data/spptreeSmithBrownShort.cr.tre")
