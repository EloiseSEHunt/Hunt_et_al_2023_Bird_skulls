### Preliminary quick ANOVAS ###
# MANOVA for comparing skull shape with trait categories
### #Adapting by EH ###
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)
library(tidyverse)
# # phylogenetic tree
sample_trees <- read.nexus("C://Users/elois/OneDrive/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C://Users/elois/OneDrive/Desktop/Ch1_analyses/EHsubsample.csv", row.names = 1)
# Habitat density factor
fact_habdens <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat)
# Diet factor
fact_diet <- species_data %>% dplyr::select(Trophic_Niche)
# Lifestyle factor
fact_life <- species_data %>% dplyr::select(Primary_Lifestyle)

# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions

load("C://Users/elois/OneDrive/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds)
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)

cs<-Y.gpa$Csize
write.csv(Y.gpa$Csize, file = "C://Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C://Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv", row.names = 1)
fact_Size <- fact_size %>% dplyr::select(x)
coords<-Y.gpa$coords
dataset<-two.d.array(coords)

# Prepare data for analyses > check that data is in the same order as tips in the tree.
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]),
            habdens=as.factor(fact_habdens[current_tree$tip.label,]),
            mig=as.factor(fact_mig[current_tree$tip.label,]),
            dev=as.factor(fact_dev[current_tree$tip.label,]),
            hab=as.factor(fact_hab[current_tree$tip.label,]),
            diet=as.factor(fact_diet[current_tree$tip.label,]),
            life=as.factor(fact_life[current_tree$tip.label,]),
            size=as.numeric(log(fact_size[current_tree$tip.label,])))


#Need to make it a dataframe of full dataset
gdf.all <- geomorph.data.frame(dataset, Csize = fact_size$x, habdens = species_data2$Habitat_Density, mig = species_data2$Migration,devmod = species_data$Dev_Mod_all, hab = species_data$Habitat, diet = species_data$Trophic_Niche, life = species_data$Primary_Lifestyle)
species_data %>%  count(Trophic_Niche, Primary_Lifestyle, sort = TRUE)
attributes(gdf.all)

habdensMVA <-procD.pgls(dataset ~ habdens, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(habdensMVA)

migMVA <-procD.pgls(dataset ~ mig, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(migMVA)

#try all together
allMVA <-procD.pgls(dataset ~ Csize*habdens*mig*diet*devmod, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(allMVA)

devmodMVA <-procD.pgls(dataset ~ devmod, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(devmodMVA)

miglifeMVA <-procD.pgls(dataset ~ Csize*devmod*hab*habdens, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(miglifeMVA)

habhabdensMVA <-procD.pgls(dataset ~ Csize*hab*devmod*life, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(habhabdensMVA)

dietMVA <-procD.pgls(dataset ~ diet*mig, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(dietMVA)

devdietMVA <-procD.pgls(dataset ~ diet*devmod, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(devdietMVA)

lifeMVA <-procD.pgls(dataset ~ life*habdens*mig*devmod, phy = current_tree, data = gdf.all, iter = 1000, SS.type = c("II"), print.progress = FALSE)
summary(lifeMVA)

### Longer mvMORPH MANOVAs ###
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)
# # phylogenetic tree
sample_trees <- read.nexus("C://Users/eloih/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C://Users/eloih/Desktop/Ch1_analyses///EHsubsample.csv", row.names = 1)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)
# Body size
fact_size <- read.csv("C://Users/eloih/Desktop/Ch1_analyses/csize.csv", row.names = 1)

# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
load("C://Users/eloih/Desktop/Ch1_analyses/slidbirds354.R")
# Prepare data for analyses > check that data is in the same order as tips in the tree.
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]),
            hab=as.factor(fact_hab[current_tree$tip.label,]),
            mig=as.factor(fact_mig[current_tree$tip.label,]),
            dev=as.factor(fact_dev[current_tree$tip.label,]),
            size=as.numeric(log(fact_size[current_tree$tip.label,])))
# fit linear model > Penalized likelihood with Pagel's lambda
# blas_set_num_threads(6) # for specifying to use multithreading through library(RhpcBLASctl)
fit <- mvgls(data~size*hab*mig*dev, data=dat, tree=current_tree, model="lambda", method = "PL")
# Multivariate test > Pillai trace
multivariate_test <- manova.gls(fit, nperm=1000, test="Pillai", nbcores=64, verbose=TRUE, type="II")
print(multivariate_test)
# save results
results <- list(fit=fit, test=multivariate_test)
#names_id <- paste("./results/tree1_",id_job,"extant_manova.Rdata", sep = "")
save(results, file="bird_manova_results.Rdata")
# end

