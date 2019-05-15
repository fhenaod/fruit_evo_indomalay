

suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(DataCombine))


# Desktop
setwd('D:/Box Sync/Projects/Projects (active)/Indo archipelago survey/Flora Malesiana analysis/MS Tree fruit evolution/Analysis')
# Laptop
#setwd('C:/Users/JB/Box Sync/Projects/Projects (active)/Indo archipelago survey/Flora Malesiana analysis/MS Tree fruit evolution/Analysis')


#--------- TRAIT DATA -----------------------------

#---- Load data and functions
rm(list = ls())
plant.data <- read.csv("spptraits0.csv")


#---- Clean the trait data
#plant.data[is.na(plant.data)] <- 0 # change blank entries to 0
plant.data$Gen_spp <- paste(plant.data$Genus, plant.data$Species, sep="_")

# Correct mis-spellings and epithets that don't play well with R
plant.data$Family <- gsub("Fagaeae", "Fagaceae", plant.data$Family)
plant.data$Gen_spp <- gsub("Ficus_d'albertisii", "Ficus_dalbertisii", plant.data$Gen_spp)

# Combine Palawan into Sundaland
plant.data$SundaPalawan <- pmax(plant.data$Sundaland, plant.data$Phil_Palawan)
plant.data$SundaPalawan <- ifelse(plant.data$Sundaland == 1, 1, plant.data$SundaPalawan)

# Remove lianas and non-woody species
plant.data <- subset(plant.data, Growth_form != "herb")
plant.data <- subset(plant.data, Growth_form != "climber")
plant.data <- subset(plant.data, Growth_form != "herb-climber")
plant.data <- subset(plant.data, Growth_form != "liana")
plant.data <- subset(plant.data, Growth_form != "aquatic-annual")

# Remove species with unknown growth form
toBeRemoved <- which(plant.data$Growth_form == "")
plant.data <- plant.data[-toBeRemoved,]
plant.data$Growth_form <- droplevels(plant.data$Growth_form)

# If there's only one length measurement, make sure it's in the "max" column
plant.data$Fruit_lg_max[is.na(plant.data$Fruit_lg_max)] <- -999
plant.data$Fruit_lg_max <- ifelse(plant.data$Fruit_lg_max == -999, plant.data$Fruit_lg_min, plant.data$Fruit_lg_max)
plant.data$Num_seeds_max[is.na(plant.data$Num_seeds_max)] <- -999
plant.data$Num_seeds_max <- ifelse(plant.data$Num_seeds_max == -999, plant.data$Num_seeds_min, plant.data$Num_seeds_max)
plant.data$Seed_lg_max[is.na(plant.data$Seed_lg_max)] <- -999
plant.data$Seed_lg_max <- ifelse(plant.data$Seed_lg_max == -999, plant.data$Seed_lg_min, plant.data$Seed_lg_max)

# Remove known highland-restricted species
plant.data$Elev_min[is.na(plant.data$Elev_min)] <- 0
plant.data <- subset(plant.data, Elev_min <= 1000)

# Remove particular taxa
plant.data <- subset(plant.data, Genus != "Gnetum") 
plant.data <- subset(plant.data, Family != "Dipterocarpaceae") 
plant.data <- subset(plant.data, Family != "Fagaceae") 

# Remove fruit types that aren't animal dispersed
plant.data <- subset(plant.data, Fruit_type != "samara") 
plant.data <- subset(plant.data, Fruit_type != "nut") 

# New column for single- vs several-seeded
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "berry", 1, NA)
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "drupe", 0, plant.data$severalseeds)
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "capsule", 1, plant.data$severalseeds)
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "pod", 1, plant.data$severalseeds)
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "nut", 0, plant.data$severalseeds)
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "syconium", 1, plant.data$severalseeds)
plant.data$severalseeds <- ifelse(plant.data$Fruit_type == "berry-dry", 1, plant.data$severalseeds)
plant.data$severalseeds <- ifelse(plant.data$Num_seeds_max > 1, 1, plant.data$severalseeds)

# Combine levels for fruit color
# see section 3.1 of Sinnott-Armstrong etal 2018 (Global Ecol Biogeogr)
plant.data$color <- ifelse(plant.data$Fruit_color == "grey-green", "goby", NA)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "red", "red", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "yellow-orange", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "yellow-red", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "black-red", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "black", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "green", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-black", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "red-purple", "red", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-red", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "yellow", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white-green", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "green-orange", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "orange", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "blue-purple", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "blue", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "black-blue", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "black-purple", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white-purple", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "grey", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-green", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "purple", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "orange-red", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-purple", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "blue-red", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white-red", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "yellow-green", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "green-red", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-yellow", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "green-purple", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "grey-brown", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "grey-blue", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-orange", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white-brown", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white-yellow", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "orange-purple", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-white", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "black-orange", "bbw", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "yellow-purple", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "orange-yellow", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "green-yellow", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "brown-olive", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "green-blue", "goby", plant.data$color)
plant.data$color <- ifelse(plant.data$Fruit_color == "white-blue", "bbw", plant.data$color)
# fruit color 1: black/blue/white
plant.data$color.bbw <- ifelse(plant.data$color == "bbw", 1, 0)
plant.data$color.bbw <- ifelse(plant.data$color == "", NA, plant.data$color.bbw)
# fruit color 2: green/orange/brown/yellow (putatively mammal-dispersed)
plant.data$color.goby <- ifelse(plant.data$color == "goby", 1, 0)
plant.data$color.goby <- ifelse(plant.data$color == "", NA, plant.data$color.goby)
# fruit color 3: red
plant.data$color.red <- ifelse(plant.data$color == "red", 1, 0)
plant.data$color.red <- ifelse(plant.data$color == "", NA, plant.data$color.red)


#---- Assemble a smaller 'working' dataset
# working dataset
spptraits <- cbind.data.frame(family=plant.data$Family, genus=plant.data$Genus, Gen_spp=plant.data$Gen_spp,
	sunda=plant.data$SundaPalawan, sulawesi=plant.data$Sulawesi, maluku=plant.data$Moluccas, 
	newguinea=plant.data$New_Guinea, fruit_lg=plant.data$Fruit_lg_max, seed_lg=plant.data$Seed_lg_max,
	num_seeds=plant.data$Num_seeds_max, severalseeds=plant.data$severalseeds,
	color_red=plant.data$color.red, color_bbw=plant.data$color.bbw, color_goby=plant.data$color.goby)
write.csv(spptraits, "spptraits1.csv", row.names=F)





#################################################################################################
#------------- PHYLOGENY ---------------------------------------------------------------------------
spptree00 <- read.tree("Zanne_Vascular_Plants_rooted.dated.tre")

ape::is.ultrametric(spptree00)
ape::is.rooted(spptree00)
ape::is.binary.tree(spptree00)

#---- Remove genera from the phylogeny that aren't in the trait database
# genera in the trait database
gentraits <- data.frame(genus=unique(spptraits$genus))
gentraits$genus <- as.character(gentraits$genus)

# genera in the phylogeny
gentree <- strsplit(spptree00$tip.label, "_")
gentree <- t(data.frame(gentree))
gentree <- data.frame(genus=unique(gentree[,1]))
gentree$genus <- as.character(gentree$genus)

# remove genera in the phylogeny that aren't in the trait database
genremove <- dplyr::anti_join(gentree, gentraits, by = "genus")
spptree0 <- spptree00
for(i in 1:nrow(genremove)){
	#i=7854
	spptree0 <- drop.tip(spptree0, tip=grep(genremove[i,1], spptree0$tip.label, value=T))	}

# check that that worked
gentree <- strsplit(spptree0$tip.label, "_")
gentree <- t(data.frame(gentree))
gentree <- data.frame(genus=unique(gentree[,1]))
gentree$genus <- as.character(gentree$genus)
genremove <- dplyr::anti_join(gentree, gentraits, by = "genus")
nrow(genremove) # should be 0

#---- Make it ultrametric (THIS RUNS OVERNIGHT!)
spptree1 <- ape::chronos(spptree0, lambda = 0) 
ape::is.ultrametric(spptree1)


write.tree(spptree, file="spptreeZanne.tree")


#---------------------------- FURTHER MANIPULATIONS --------------------------------------------
#rm(list = ls())
#spptree1 <- spptree
spptree1 <- ape::read.tree("spptreeZanne.tree") 
spptraits1 <- read.csv("spptraits1.csv")


#---- Remove spp from phylogeny that ARE in the traits database, but keep their genera 
# species in the trait database
speciestraits <- data.frame(species=unique(spptraits1$Gen_spp))
speciestraits$species <- as.character(speciestraits$species)

# Rename spp in the phylogeny that are in the traits database
tips <- data.frame(species=spptree1$tip.label, stringsAsFactors=F)
names(tips) <- "species"
tips$species <- as.character(tips$species)

# species in the phylogeny that are in the traits database
sppboth <- dplyr::inner_join(speciestraits, tips, by="species")

# Rename them and insert them back into the phylogeny
sppboth$namesnew <- paste(sppboth$species, "1", sep="")
tipsnew <- DataCombine::FindReplace(data=tips, Var="species", replaceData=sppboth, 
	from="species", to="namesnew", exact=TRUE, vector=TRUE)
spptree1$tip.label <- tipsnew


#---- Add species in the traits database to the phylogeny
# species in the trait database
SppToAdd <- data.frame(species=unique(spptraits1$Gen_spp))
SppToAdd$species <- as.character(speciestraits$species)

# Can't add species that don't have corresponding genera in the phylogeny
# genera in the phylogeny
gentree <- strsplit(spptree1$tip.label, "_")
gentree <- t(data.frame(gentree))
gentree <- data.frame(genus=unique(gentree[,1]))
gentree$genus <- as.character(gentree$genus) 

# adding a 'genus' column to SppToAdd 
tmp <- strsplit(SppToAdd$species, "_")
tmp <- t(data.frame(tmp))
tmp <- data.frame(genus=tmp[,1])
tmp$genus <- as.character(tmp$genus)
SppToAdd$genus <- tmp$genus

# remove species from SppToAdd that don't have a genus in the phylogeny
SppToAdd <- dplyr::inner_join(SppToAdd, gentree, by="genus")


# Add species to the phylogeny (at the root node of the appropriate genus) that are in the traits database but not the phylogeny
spptree <- spptree1
st <- Sys.time()
for(i in 1:nrow(SppToAdd)){
   tryCatch({
	#i=169
	spptree <- phytools::add.species.to.genus(spptree, SppToAdd[i,1], where="root")
	#ape::is.ultrametric(spptree)
   }, error=function(e){}) # end of the tryCatch function
}
et <- Sys.time()
et-st
ape::is.ultrametric(spptree)


# Remove species in phylogeny that are not in trait data
in.spec <- spptree$tip.label %in% spptraits1$Gen_spp
spptree <- drop.tip(spptree, spptree$tip.label[!in.spec])

unique(spptree$node.label)

write.tree(spptree, file="spptreeZanneShort.tree")



















#################################################################################################
#------------ MORE PHYLOGENY MANIPULATIONS - NOT CURRENTLY USED --------------------------------------
# Deal with any negative branch lengths
summary(spptree00$edge.length) # There's one branch with 'NaN' branch length
edges <- spptree00$edge.length
edges[is.nan(edges)] <- 0 # replace NaN values with 0 
#edges[edges == 0] <- 0.0001 # replace the 0 branch lengths with very small, but positive, branch lengths
spptree0 <- ape::compute.brlen(spptree00, edges) # replace the branch lengths in the phylogeny with a vector of corrected lengths

# Make the tree dichotomous
spptree0 <- ape::multi2di(spptree0, random = TRUE) # resolves polytomies by randomly assigning dichotomous breakpoints
spptree0 <- ape::collapse.singles(spptree0, root.edge=TRUE) # resolve nodes with only 1 (instead of 2) descendents
ape::is.binary.tree(spptree0)

# Provide a time calibration (age at which outgroup, Amborella, split from the rest)
mycalibration <- makeChronosCalib(spptree0, node="root", age.max=160) # age.max in millions of years

# Loop to find optimum lambda value from Cadotte & Davies ("Phylogenies in Ecology", p.35)
l <- c(0:10, seq(from=12, to=100, by=2))
LL.out <- NULL
Sys.time()
for(i in 1:length(l)){LL.out[i] <- attributes(ape::chronos(spptree0, lambda=l[i], model="relaxed", calibration=mycalibration, control=chronos.control()))$ploglik	}
Sys.time() # took 36 hours on my laptop
l[LL.out==max(LL.out)] # find the lambda value for which log likelihood is maximized
plot(log(l+1), LL.out, type="l", lwd=3, col="gray", xlab=expression(paste("Ln(",lambda,")")), ylab="Log likelihood")

# Save lambda optimization results
output <- cbind.data.frame(l, LL.out)
names(output) <- c("lambda", "log.likelihood")
write.csv(output, "lambda_likelihoods.csv", row.names=F)
#-------------------------------------------------------------------------------------------------------



















############ NOT CURRENTLY USED (I DID THIS IN ACCESS INSTEAD) ##########################################
# Combine duplicate rows
# Note: This removes the character & factor columns
plant.data1 <- plyr::ddply(plant.data, "Gen_spp", numcolwise(max))























