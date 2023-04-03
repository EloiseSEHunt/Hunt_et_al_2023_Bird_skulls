### Plotting traits on tree_Figure 1 ###
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)
library(tidyverse)
library(evobiR)

# Phylogenetic tree
sample_trees <- read.nexus("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
load("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/slidbirds354.R")

Y.gpa<-gpagen(slidbirds) #generalised procrustes analysis
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
# Ecological factors
species_data <- read.csv("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/EHsubsample.csv")
row.names(species_data) <- species_data$Tip_Label

species_data<-ReorderData(current_tree, species_data, taxa.names="Tip_Label")

#Make a color palatte for each trait (with 3 characters in each trait)
Hab_col <- c("#8B6508","#FFD700","#FFF8DC")
Mig_col <- c("#E69F00", "#56B4E9", "#009E73")
Dev_col <- c("#0072B2", "#BF3EFF", "#FF1493")

# Sets up more space to right of figure for legend
par(mar=c(5, 4, 4, 8), xpd=TRUE)

# Plot the fan tree with offset tip labels
plot(current_tree,type="fan",show.tip.label=TRUE, label.offset = 5.5, cex = 0.3)

#Add the traits
Habdens<-as.factor(species_data$Habitat_Density)
Habdens2<-to.matrix(Habdens,levels(Habdens))
rownames(Habdens2)<-species_data$Tip_Label
Habdens2<-Habdens2[current_tree$tip.label,]
tiplabels(pie=Habdens2,piecol=Hab_col,cex=0.2,offset = 1)

Mig<-as.factor(species_data$Migration)
Mig2<-to.matrix(Mig,levels(Mig))
rownames(Mig2)<-species_data$Tip_Label
Mig2<-Mig2[current_tree$tip.label,]
tiplabels(pie=Mig2,piecol=Mig_col,cex=0.2,offset = 2.5)

Devmod<-as.factor(species_data$Dev_Mod_all)
Devmod2<-to.matrix(Devmod,levels(Devmod))
rownames(Devmod2)<-species_data$Tip_Label
Devmod2<-Devmod2[current_tree$tip.label,]
tiplabels(pie=Devmod2,piecol=Dev_col,cex=0.2,offset = 4)

# Legend
legend("topright",inset=c(-0.3, 0), levels(Habdens), title="Habitat density", pch=21,pt.bg=Hab_col, pt.cex=1, cex=0.5, box.col="white")
legend("right", inset=c(-0.3, 0), levels(Mig),title="Migration", pch=21,pt.bg=Mig_col, pt.cex=1, cex=0.5, box.col="white")
legend("bottomright",inset=c(-0.3, 0),levels(Devmod),title="Developmental mode", pch=21,pt.bg=Dev_col, pt.cex=1, cex=0.5, box.col="white")

### Plotting morphospaces_Figure 2 ###
library(mvMORPH)
library(parallel)
library(dplyr)
library(geiger)
library(geomorph)

# # phylogenetic tree
sample_trees <- read.nexus("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/correct.composite.phylogeny.oct.26.nex")
# Ecological factors
species_data <- read.csv("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/EHsubsample.csv", row.names = 1)
# Habitat factor
fact_hab <- species_data %>% dplyr::select(Habitat_Density)
# Migration factor
fact_mig <- species_data %>% dplyr::select(Migration)
# Developmental mode factor
fact_dev <- species_data %>% dplyr::select(Dev_Mod_all)

# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
load("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/slidbirds354.R")
Y.gpa<-gpagen(slidbirds) #generalised procrustes analysis
nc1<-name.check(sample_trees,two.d.array(slidbirds))
nc1$data_not_tree
current_tree<-drop.tip(sample_trees, nc1$tree_not_data)
write.csv(Y.gpa$Csize, file = "C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/csize.csv")
fact_size <- read.csv("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/csize.csv", row.names = 1)
proOut<-Y.gpa$coords
dataset<-two.d.array(coords)
str(coords)
PCA.3D <- geomorph::gm.prcomp( proOut )

# Basic ggplot PCA with convex hulls
library( ggplot2 )
library (ggConvexHull)
library (ggpubr)
library( geomorph )

#read in csv for categories for Habitat Density
classifier <- read.csv("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/EHsubsample.csv")
str(classifier)
classifier$Habitat_density
# import family information from classifier csv file
Habitat_density <- factor( classifier$Habitat_density ) # create factor
lev.group <- length( levels( Habitat_density ) )
Habitat_density <-as.factor(Habitat_density)
# set colours, one for each group 
col.gp.1 <- c("darkorange3", "darkcyan",  "darkgoldenrod1")
names( col.gp.1 )   <- levels( Habitat_density )
# Assign one colour to each specimen
col.gp   <- col.gp.1[ match( Habitat_density, names( col.gp.1 ) ) ]

# Set PCA axis labels
xlab <- paste("Principal Component 1 (",signif((PCA.3D$d/sum(PCA.3D$d)*100),3)[1] , "%)", sep = "")
ylab <- paste("Principal Component 2 (",signif((PCA.3D$d/sum(PCA.3D$d)*100),3)[2] , "%)", sep = "")

# create ggplot dataframe
# classifier is a csv file with one row for each specimen, and one column for each category
PC1 <-data.frame(cbind(PCA.3D$x[,c(1:14)],classifier,col.gp,dimnames(proOut)[[3]]))
length(col.gp)
# Plot PCs 1 and 2 for PCA
P1 <-ggplot(PC1, aes(x=Comp1 ,y=Comp2,colour = Habitat_density, fill = Habitat_density) )+
  geom_convexhull(alpha = 0.2, show.legend = TRUE)+
  geom_point(colour = "black", shape=21, size=3.5,show.legend = TRUE)+
  labs(x = xlab, y= ylab)+
  scale_fill_manual(values = col.gp.1)+
  scale_colour_manual(values = col.gp.1)+
  theme_bw()
plot(P1)

p1a <- P1 + theme(text = element_text(size = 15))
plot(p1a)

#read in csv for categories for Migration
classifier <- read.csv("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/EHsubsample.csv")
str(classifier)
classifier$Migration
# import family information from classifier csv file
Migration <- factor( classifier$Migration ) # create factor
lev.group <- length( levels( Migration ) )
Migration <-as.factor(Migration)
# set colours, one for each group 
col.gp.1 <- c("darkorange2", "orchid1",  "purple3") 
names( col.gp.1 )   <- levels( Migration )
# Assign one colour to each specimen
col.gp   <- col.gp.1[ match( Migration, names( col.gp.1 ) ) ]

# Set PCA axis labels
xlab <- paste("Principal Component 1 (",signif((PCA.3D$d/sum(PCA.3D$d)*100),3)[1] , "%)", sep = "")
ylab <- paste("Principal Component 2 (",signif((PCA.3D$d/sum(PCA.3D$d)*100),3)[2] , "%)", sep = "")

# create ggplot dataframe
# classifier is a csv file with one row for each specimen, and one column for each category
PC1 <-data.frame(cbind(PCA.3D$x[,c(1:14)],classifier,col.gp,dimnames(proOut)[[3]]))
length(col.gp)
# Plot PCs 1 and 2 for PCA
P2 <-ggplot(PC1, aes(x=Comp1 ,y=Comp2,colour = Migration, fill = Migration) )+
  geom_convexhull(alpha = 0.2, show.legend = TRUE)+
  geom_point(colour = "black", shape=21, size=3.5,show.legend = TRUE)+
  labs(x = xlab, y= ylab)+
  scale_fill_manual(values = col.gp.1)+
  scale_colour_manual(values = col.gp.1)+
  theme_bw()
plot(P2)

P2a <- P2 + theme(text = element_text(size = 15))
plot(P2a)

#read in csv for categories for Dev Mod
classifier <- read.csv("C:/Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses/EHsubsample.csv")
str(classifier)
classifier$Dev_Mod
# import family information from classifier csv file
Dev_Mod <- factor( classifier$Dev_Mod ) # create factor
lev.group <- length( levels( Dev_Mod ) )
Dev_Mod <-as.factor(Dev_Mod)
# set colours, one for each group 
col.gp.1 <- c("darkolivegreen3", "brown",  "skyblue3") 
names( col.gp.1 )   <- levels( Dev_Mod )
# Assign one colour to each specimen
col.gp   <- col.gp.1[ match( Dev_Mod, names( col.gp.1 ) ) ]

# Set PCA axis labels
xlab <- paste("Principal Component 1 (",signif((PCA.3D$d/sum(PCA.3D$d)*100),3)[1] , "%)", sep = "")
ylab <- paste("Principal Component 2 (",signif((PCA.3D$d/sum(PCA.3D$d)*100),3)[2] , "%)", sep = "")

# create ggplot dataframe
# classifier is a csv file with one row for each specimen, and one column for each category
PC1 <-data.frame(cbind(PCA.3D$x[,c(1:14)],classifier,col.gp,dimnames(proOut)[[3]]))
length(col.gp)
# Plot PCs 1 and 2 for PCA
P3 <-ggplot(PC1, aes(x=Comp1 ,y=Comp2,colour = Dev_Mod, fill = Dev_Mod) )+
  geom_convexhull(alpha = 0.2, show.legend = TRUE)+
  geom_point(colour = "black", shape=21, size=3.5,show.legend = TRUE)+
  labs(x = xlab, y= ylab)+
  scale_fill_manual(values = col.gp.1)+
  scale_colour_manual(values = col.gp.1)+
  theme_bw()
plot(P3)

P3a <- P3 + theme(text = element_text(size = 15))
plot(P3a)

# Plot histogram of proportion of PC variances
summary(PCA.3D)
par(bg = 'white')
pvar <- (PCA.3D$sdev^2)/(sum(PCA.3D$sdev^2)) ###graph of PC proportional variance
names(pvar) <- seq(1:length(pvar))
barplot(pvar, xlab= "Principal Components", ylab = "Proportion of Variance", ylim = c(0,1))

### Plotting rates_Figure 3 ###
library(dplyr)
library(ggplot2)
library(tidyverse)

# CSV with rates
comp <- read.csv("C://Users/elois/OneDrive - Natural History Museum/Desktop/Ch1_analyses///rates.csv")

p <- ggplot(comp, aes(x = factor(Category, level = c('Dense', 'Semi-open', 'Open', 'Non-migratory', 'Partially migratory', 'Migratory', 'Precocial', 'Semi-precocial', 'Altricial')), y = Rates)) +
  geom_point(aes(color = Factor), size = 10) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_bw() +
  theme(
    legend.position = "bottom" ,
    legend.direction = "horizontal",
    axis.text.x = element_text (angle = 45, hjust = 1) ,
    axis.title.x = element_blank() ,
    text = element_text(size = 20) ,
    # Hide panel borders and remove grid lines
    panel.border = element_blank()) +
  labs(y = "Evolutionary rate")

print (p)