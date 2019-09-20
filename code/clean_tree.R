library(ape)
library(castor)
spptree1<-read.tree("data/spptreeZanneShort.tree") 

is.binary(spptree1)
spptree2<-multi2di(spptree1)
min(spptree2$edge.length)
spptree2$edge.length[spptree2$edge.length<=0] <-1e-6

is.ultrametric(spptree2)
zanne_tree<-castor:::extend_tree_to_height(spptree2)$tree
is.ultrametric(zanne_tree)
is.binary(zanne_tree)

write.tree(zanne_tree, file = "data/zanne_tree.tre")
