###Overall skull rates###

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
# Developmental mode (no family level inferences) factor
fact_devnofam <- species_data %>% dplyr::select(Dev_Mod_no_family_order)


load("C://Users/eloih/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds)
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
Y.gpa<-gpagen(slidbirds)
#cs<-Y.gpa$Csize
write.csv(Y.gpa$Csize, file = "C://Users/eloih/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C://Users/eloih/Desktop/Ch1_analyses/csize.csv", row.names = 1)
coords<-Y.gpa$coords
dataset<-two.d.array(coords)

########## habitat density rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_hab[current_tree$tip.label,])
names(cat_eco) = current_tree$tip.label
# i)
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?
# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
#Add the function to reconstruct from a sample of trees 
paintAllTree <- function(tree, ancestral, tips){  
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }  
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }  
  return(treebis)
}

# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. 
#Ideally this should be performed on a sample of stochastic maps.
#Or by using a mapping constructed from the ML reconstruction for instance
simm_tr<- make.simmap(current_tree, model="ARD", cat_eco, nsim=100)
a_single_simmap_tree <- paintAllTree(current_tree, describe.simmap(simm_tr), as.character(cat_eco))


dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=a_single_simmap_tree, model="BMM", method = "PL", error=TRUE)
# if you want to add a covariate, e.g. size, then use the following instead:
# fit <- mvgls(data~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)
# Save the results.
# Multiple options here again. If you just want the rates (and not say, the evolutionary covariances), you can save a lot of memory space by saving only that:
#name_rates_file <- paste("estimated_rates_tree_",id_job,".txt",sep="")
#write.table(fit$param, file=name_rates_file)
# # or save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object!
results <- list(fit=fit$param, tree=a_single_simmap_tree)
#names_id <- paste("tree1_75_80_",id_job,"diet_size_BMM_.Rdata", sep = "")
save(results, file="bird_hab_rates.Rdata")


########## migration rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_mig[current_tree$tip.label,])
names(cat_eco) = current_tree$tip.label
# i)
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?
# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
#Add the function to reconstruct from a sample of trees 
paintAllTree <- function(tree, ancestral, tips){  
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }  
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }  
  return(treebis)
}

# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. 
#Ideally this should be performed on a sample of stochastic maps.
#Or by using a mapping constructed from the ML reconstruction for instance
simm_tr<- make.simmap(current_tree, model="ARD", cat_eco, nsim=100)
a_single_simmap_tree <- paintAllTree(current_tree, describe.simmap(simm_tr), as.character(cat_eco))


dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=a_single_simmap_tree, model="BMM", method = "PL", error=TRUE)
# if you want to add a covariate, e.g. size, then use the following instead:
# fit <- mvgls(data~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)
# Save the results.
# Multiple options here again. If you just want the rates (and not say, the evolutionary covariances), you can save a lot of memory space by saving only that:
#name_rates_file <- paste("estimated_rates_tree_",id_job,".txt",sep="")
#write.table(fit$param, file=name_rates_file)
# # or save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object!
results <- list(fit=fit$param, tree=a_single_simmap_tree)
#names_id <- paste("tree1_75_80_",id_job,"diet_size_BMM_.Rdata", sep = "")
save(results, file="bird_mig_rates.Rdata")
results$fit

########## developmental mode rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_dev[current_tree$tip.label,])
names(cat_eco) = current_tree$tip.label
# i)
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?
# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
#Add the function to reconstruct from a sample of trees 
paintAllTree <- function(tree, ancestral, tips){  
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }  
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }  
  return(treebis)
}

# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. 
#Ideally this should be performed on a sample of stochastic maps.
#Or by using a mapping constructed from the ML reconstruction for instance
simm_tr<- make.simmap(current_tree, model="ARD", cat_eco, nsim=100)
a_single_simmap_tree <- paintAllTree(current_tree, describe.simmap(simm_tr), as.character(cat_eco))


dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=a_single_simmap_tree, model="BMM", method = "PL", error=TRUE)
# if you want to add a covariate, e.g. size, then use the following instead:
# fit <- mvgls(data~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)
# Save the results.
# Multiple options here again. If you just want the rates (and not say, the evolutionary covariances), you can save a lot of memory space by saving only that:
#name_rates_file <- paste("estimated_rates_tree_",id_job,".txt",sep="")
#write.table(fit$param, file=name_rates_file)
# # or save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object!
results <- list(fit=fit$param, tree=a_single_simmap_tree)
#names_id <- paste("tree1_75_80_",id_job,"diet_size_BMM_.Rdata", sep = "")
save(results, file="bird_dev_rates.Rdata")



### SES ###
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)

# # phylogenetic tree
sample_trees <- read.nexus("C://Users/elois/OneDrive/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C://Users/elois/OneDrive/Desktop/Ch1_analyses///EHsubsample.csv", row.names = 1)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)
# Developmental mode (no family level inferences) factor
fact_devnofam <- species_data %>% dplyr::select(Dev_Mod_no_family_order)


load("C://Users/elois/OneDrive/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds)
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
Y.gpa<-gpagen(slidbirds)
#cs<-Y.gpa$Csize
write.csv(Y.gpa$Csize, file = "C://Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C://Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv", row.names = 1)
coords<-Y.gpa$coords
dataset<-two.d.array(coords)

########## habitat density rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_hab[current_tree$tip.label,])
names(cat_eco) = current_tree$tip.label
# i)
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?
# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
#Add the function to reconstruct from a sample of trees 
paintAllTree <- function(tree, ancestral, tips){  
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }  
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }  
  return(treebis)
}

# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. 
#Ideally this should be performed on a sample of stochastic maps.
#Or by using a mapping constructed from the ML reconstruction for instance
simm_tr<- make.simmap(current_tree, model="ARD", cat_eco, nsim=100)
a_single_simmap_tree <- paintAllTree(current_tree, describe.simmap(simm_tr), as.character(cat_eco))

# get the effect size (SES)
ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
ls()
results$fit

# get table with all results
pillai <- results$test$stat
terms <- results$test$terms
p_val <- results$test$pvalue
ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
datmat<-rbind(pillai,ses,p_val)
colnames(datmat)<-terms
data_tmp <- as_tibble(datmat)
data_tmp <- data_tmp %>% mutate(result = c("Pallai's Test Statistic","SES", "P value"))
data <- bind_rows(data,data_tmp) 
ls()

### Module rates ###

library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)

# # phylogenetic tree
sample_trees <- read.nexus("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/EHsubsample.csv", row.names = 1)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds) #generalised procrustes analysis
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
write.csv(Y.gpa$Csize, file = "C:/Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv", row.names = 1)
proOut<-Y.gpa$coords
#dataset<-two.d.array(coords)


#####Start new code for module rates ######

#Whole Skull
species_data<- species_data[dimnames(proOut)[[3]], ]
fulldf <- geomorph.data.frame(shape = proOut, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.lm(shape ~ habdens, data = fulldf)

#make moduledef matrix
Moduledefs <-  read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/module.definitions.aug.2022.csv")
left.lm<-c(23,24,25,26,27,28,29,30,31,32,33,36,34)
Moduledefs <- Moduledefs[-left.lm,]

#each module
rostrumcoords <- proOut[which(Moduledefs$posthoc == 1), , ]
rostrumdf <- geomorph.data.frame(shape = rostrumcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = rostrumdf)

vaultcoords <- proOut[which(Moduledefs$posthoc == 2), , ]
vaultdf <- geomorph.data.frame(shape = vaultcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = vaultdf)

sphenoidcoords <- proOut[which(Moduledefs$posthoc == 3), , ]
sphenoiddf <- geomorph.data.frame(shape = sphenoidcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = sphenoiddf)

palatecoords <- proOut[which(Moduledefs$posthoc == 4), , ]
palatedf <- geomorph.data.frame(shape = palatecoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = palatedf)

jointcoords <- proOut[which(Moduledefs$posthoc == 5), , ]
jointdf <- geomorph.data.frame(shape = jointcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = jointdf)

nariscoords <- proOut[which(Moduledefs$posthoc == 6), , ]
narisdf <- geomorph.data.frame(shape = nariscoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = narisdf)

occipitalcoords <- proOut[which(Moduledefs$posthoc == 7), , ]
occipitaldf <- geomorph.data.frame(shape = occipitalcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ habdens, phy = current_tree, data = occipitaldf)

## developmental mode rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_hab[current_tree$tip.label,]) #was dev_mig
names(cat_eco) = current_tree$tip.label
# i)
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?
# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
# paintAllTree <- function(tree, ancestral, tips){  
#   
#   if(inherits(ancestral, "describe.simmap")){
#     names_grps <- colnames(ancestral$ace)
#     statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
#   }else{
#     names_grps <- colnames(ancestral$lik.anc)
#     statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
#   }  
#   
#   combined = as.character(c(tips, statesNodes))
#   treebis=tree
#   for(i in sort(tree$edge[,2])){
#     treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
#   }  
#   return(treebis)
# }
# 
# # First estimate ancestral state using ML in "ape"
# ace_habitat <- ace(cat_eco, sample_trees[[id_job]], type = "discrete", model="ARD")
# # convert the ML reconstruction to a "simmap" tree object
# tree_simmap <- paintAllTree(sample_trees[[id_job]], ace_habitat, as.character(cat_eco))
# # iii) marginal reconstructions by stochastic mapping
# # i.e. do the same using stochastic mapping instead of ML
# nsim = 100 # Number of stochastic mappings
# my_trees<-make.simmap(sample_trees[[id_job]], cat_eco , model="ARD", nsim=nsim)
# 
# # Use the collection of stochastic mappings to reconstruct the marginal ancestral states
# distrib <-summary(my_trees, plot=TRUE)
# # Use these reconstruction to build a simmap tree
# tree_simmap  <- paintAllTree(sample_trees[[id_job]], distrib, as.character(cat_eco))

dataset<-two.d.array(rostrumcoords)

## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##
## !! check that the object names match the ones you used above
# Prepare data for analyses in a list
#dat <- list(data=as.matrix(dataset[sample_trees[[id_job]]$tip.label,]))
# If you want to account for size, rather than removing size first, you can use the following instead
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
# if you want to add a covariate, e.g. size, then use the following instead:
# fit <- mvgls(data~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)
# Save the results.
# Multiple options here again. If you just want the rates (and not say, the evolutionary covariances), you can save a lot of memory space by saving only that:
#name_rates_file <- paste("estimated_rates_tree_",id_job,".txt",sep="")
#write.table(fit$param, file=name_rates_file)
# # or save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object!
results <- list(fit=fit$param, tree=tree_simmap)
#names_id <- paste("tree1_75_80_",id_job,"diet_size_BMM_.Rdata", sep = "")
save(results, file="rostrum_hab_rates_new module.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/rostrum_hab_rates_new module.Rdata")
ls()
results$fit

dataset<-two.d.array(vaultcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="vault_hab_rates_new module.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/vault_hab_rates_new module.Rdata")
results$fit

dataset<-two.d.array(sphenoidcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="sphenoid_hab_rates_new module.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/sphenoid_hab_rates_new module.Rdata")
results$fit

dataset<-two.d.array(palatecoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="palate_hab_rates_new module.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/palate_hab_rates_new module.Rdata")
results$fit

dataset<-two.d.array(jointcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="joint_hab_rates.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/joint_hab_rates_new module.Rdata")
results$fit

dataset<-two.d.array(nariscoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="naris_hab_rates_new module.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/naris_hab_rates_new module.Rdata")
results$fit

dataset<-two.d.array(occipitalcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="occipital_hab_rates_new module.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/occipital_hab_rates_new module.Rdata")
results$fit

#####
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)

# # phylogenetic tree
sample_trees <- read.nexus("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/EHsubsample.csv", row.names = 1)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds) #generalised procrustes analysis
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
write.csv(Y.gpa$Csize, file = "C:/Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv", row.names = 1)
proOut<-Y.gpa$coords
#dataset<-two.d.array(coords)


#####Start new code for module rates ######

#Whole Skull
species_data<- species_data[dimnames(proOut)[[3]], ]
fulldf <- geomorph.data.frame(shape = proOut, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.lm(shape ~ migration, data = fulldf)

#make moduledef matrix
Moduledefs <-  read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/module.definitions.aug.2022.csv")
left.lm<-c(23,24,25,26,27,28,29,30,31,32,33,36,34)
Moduledefs <- Moduledefs[-left.lm,]

#each module
rostrumcoords <- proOut[which(Moduledefs$posthoc == 1), , ]
rostrumdf <- geomorph.data.frame(shape = rostrumcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = rostrumdf)

vaultcoords <- proOut[which(Moduledefs$posthoc == 2), , ]
vaultdf <- geomorph.data.frame(shape = vaultcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = vaultdf)

sphenoidcoords <- proOut[which(Moduledefs$posthoc == 3), , ]
sphenoiddf <- geomorph.data.frame(shape = sphenoidcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = sphenoiddf)

palatecoords <- proOut[which(Moduledefs$posthoc == 4), , ]
palatedf <- geomorph.data.frame(shape = palatecoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = palatedf)

jointcoords <- proOut[which(Moduledefs$posthoc == 5), , ]
jointdf <- geomorph.data.frame(shape = jointcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = jointdf)

nariscoords <- proOut[which(Moduledefs$posthoc == 6), , ]
narisdf <- geomorph.data.frame(shape = nariscoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = narisdf)

occipitalcoords <- proOut[which(Moduledefs$posthoc == 7), , ]
occipitaldf <- geomorph.data.frame(shape = occipitalcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ migration, phy = current_tree, data = occipitaldf)

#### OLD RATES METHOD ####
########## developmental mode rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_mig[current_tree$tip.label,]) #was dev_mig
names(cat_eco) = current_tree$tip.label
# i)
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?
# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
# paintAllTree <- function(tree, ancestral, tips){  
#   
#   if(inherits(ancestral, "describe.simmap")){
#     names_grps <- colnames(ancestral$ace)
#     statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
#   }else{
#     names_grps <- colnames(ancestral$lik.anc)
#     statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
#   }  
#   
#   combined = as.character(c(tips, statesNodes))
#   treebis=tree
#   for(i in sort(tree$edge[,2])){
#     treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
#   }  
#   return(treebis)
# }
# 
# # First estimate ancestral state using ML in "ape"
# ace_habitat <- ace(cat_eco, sample_trees[[id_job]], type = "discrete", model="ARD")
# # convert the ML reconstruction to a "simmap" tree object
# tree_simmap <- paintAllTree(sample_trees[[id_job]], ace_habitat, as.character(cat_eco))
# # iii) marginal reconstructions by stochastic mapping
# # i.e. do the same using stochastic mapping instead of ML
# nsim = 100 # Number of stochastic mappings
# my_trees<-make.simmap(sample_trees[[id_job]], cat_eco , model="ARD", nsim=nsim)
# 
# # Use the collection of stochastic mappings to reconstruct the marginal ancestral states
# distrib <-summary(my_trees, plot=TRUE)
# # Use these reconstruction to build a simmap tree
# tree_simmap  <- paintAllTree(sample_trees[[id_job]], distrib, as.character(cat_eco))

dataset<-two.d.array(rostrumcoords)

## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##
## !! check that the object names match the ones you used above
# Prepare data for analyses in a list
#dat <- list(data=as.matrix(dataset[sample_trees[[id_job]]$tip.label,]))
# If you want to account for size, rather than removing size first, you can use the following instead
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
# if you want to add a covariate, e.g. size, then use the following instead:
# fit <- mvgls(data~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)
# Save the results.
# Multiple options here again. If you just want the rates (and not say, the evolutionary covariances), you can save a lot of memory space by saving only that:
#name_rates_file <- paste("estimated_rates_tree_",id_job,".txt",sep="")
#write.table(fit$param, file=name_rates_file)
# # or save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object!
results <- list(fit=fit$param, tree=tree_simmap)
#names_id <- paste("tree1_75_80_",id_job,"diet_size_BMM_.Rdata", sep = "")
save(results, file="rostrum_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/rostrum_mig_rates_correct_modules.Rdata")
ls()
results$fit

dataset<-two.d.array(vaultcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="vault_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/vault_mig_rates_correct_modules.Rdata")
results$fit

dataset<-two.d.array(sphenoidcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="sphenoid_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/sphenoid_mig_rates_correct_modules.Rdata")
results$fit

dataset<-two.d.array(palatecoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="palate_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/palate_mig_rates_correct_modules.Rdata")
results$fit

dataset<-two.d.array(jointcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="joint_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/joint_mig_rates_correct_modules.Rdata")
results$fit

dataset<-two.d.array(nariscoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="naris_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/naris_mig_rates_correct_modules.Rdata")
results$fit

dataset<-two.d.array(occipitalcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="occipital_mig_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/occipital_mig_rates_correct_modules.Rdata")
results$fit

###
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)

# # phylogenetic tree
sample_trees <- read.nexus("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/EHsubsample.csv", row.names = 1)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds) #generalised procrustes analysis
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
write.csv(Y.gpa$Csize, file = "C:/Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/csize.csv", row.names = 1)
proOut<-Y.gpa$coords
#dataset<-two.d.array(coords)


#####Start new code for module rates ######

#Whole Skull
species_data<- species_data[dimnames(proOut)[[3]], ]
fulldf <- geomorph.data.frame(shape = proOut, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.lm(shape ~ devmod, data = fulldf)

#make moduledef matrix
Moduledefs <-  read.csv("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/module.definitions.aug.2022.csv")
left.lm<-c(23,24,25,26,27,28,29,30,31,32,33,36,34)
Moduledefs <- Moduledefs[-left.lm,]

#each module
rostrumcoords <- proOut[which(Moduledefs$posthoc == 1), , ]
rostrumdf <- geomorph.data.frame(shape = rostrumcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = rostrumdf)

vaultcoords <- proOut[which(Moduledefs$posthoc == 2), , ]
vaultdf <- geomorph.data.frame(shape = vaultcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = vaultdf)

sphenoidcoords <- proOut[which(Moduledefs$posthoc == 3), , ]
sphenoiddf <- geomorph.data.frame(shape = sphenoidcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = sphenoiddf)

palatecoords <- proOut[which(Moduledefs$posthoc == 4), , ]
palatedf <- geomorph.data.frame(shape = palatecoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = palatedf)

jointcoords <- proOut[which(Moduledefs$posthoc == 5), , ]
jointdf <- geomorph.data.frame(shape = jointcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = jointdf)

nariscoords <- proOut[which(Moduledefs$posthoc == 6), , ]
narisdf <- geomorph.data.frame(shape = nariscoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = narisdf)

occipitalcoords <- proOut[which(Moduledefs$posthoc == 7), , ]
occipitaldf <- geomorph.data.frame(shape = occipitalcoords, habdens = species_data$Habitat_Density, migration = species_data$Migration, devmod = species_data$Dev_Mod_all)
procD.pgls(shape ~ devmod, phy = current_tree, data = occipitaldf)


#### OLD RATES METHOD ####
########## developmental mode rates
## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##
# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)
cat_eco = as.factor(fact_dev[current_tree$tip.label,]) #was dev_mig
names(cat_eco) = current_tree$tip.label
tree_simmap <- make.simmap(current_tree, cat_eco , model="ARD", nsim=1) # replace ARD by SYM?

## rostrum
dataset<-two.d.array(rostrumcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))

fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="module_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/module_dev_rates_correct_modules.Rdata")
ls()
results$fit

## vault
dataset<-two.d.array(vaultcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="vault_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/vault_dev_rates_correct_modules.Rdata")
ls()
results$fit

## sphenoid
dataset<-two.d.array(sphenoidcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="sphenoid_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/sphenoid_dev_rates_correct_modules.Rdata")
results$fit

## palate
dataset<-two.d.array(palatecoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="palate_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/palate_dev_rates_correct_modules.Rdata")
results$fit

## joint
dataset<-two.d.array(jointcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="joint_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/joint_dev_rates_correct_modules.Rdata")
results$fit

## naris
dataset<-two.d.array(nariscoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="naris_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/naris_dev_rates_correct_modules.Rdata")
results$fit

## occipital
dataset<-two.d.array(occipitalcoords)
dat <- list(data=as.matrix(dataset[current_tree$tip.label,]), 
            size=as.numeric(fact_size[current_tree$tip.label,]))
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="occipital_dev_rates_correct_modules.Rdata")
load("C:/Users/elois/OneDrive/Desktop/Ch1_analyses/occipital_dev_rates_correct_modules.Rdata")
results$fit


