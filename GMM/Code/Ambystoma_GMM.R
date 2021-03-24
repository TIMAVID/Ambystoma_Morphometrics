## Load in data ##

download.file(   "https://github.com/TIMAVID/Ambystoma/blob/fe4c0fe3f83265c6aa656946b696096524cdd26d/GMM/Data/GMM_data.RData?raw=true",   "GMM_data.RData")
load("GMM_data.RData")

library(geomorph)


# Generalized procrustes analysis

GPA_landmarks_sub <- gpagen(GMM_data_sub$land)

#Create geomorph data frame

Amb_gdf<-geomorph.data.frame(coords=GPA_landmarks_sub$coords,
                             size=GPA_landmarks_sub$Csize, species=GMM_data_sub$species)
## PCA ##

Amb_PCA <- gm.prcomp(GPA_landmarks_sub$coords)

?gm.prcomp

### PCA vizualization ###

PC_scores <- as.data.frame(Amb_PCA$x)
library(ggplot2)
library(grid)
library(gridExtra)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
percentage <- round(Amb_PCA$sdev / sum(Amb_PCA$sdev) * 100, 2)
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

# library(RColorBrewer)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# my.colors=rainbow(28) #set up color palette of rainbow colors with n = 14
# plot(1:28, pch=19, cex=2, col=my.colors)

library(randomcoloR)
n <- 14
palette <- distinctColorPalette(n)

p<-ggplot(PC_scores,aes(x=Comp1,y=Comp2,color=GMM_data_sub$species ))
p<-p+geom_point(size =5)+theme + xlab(percentage[1]) + ylab(percentage[2]) +scale_color_manual(values = palette)
# p + stat_ellipse()
p

## ANOVA ##

# Without size
Amb_anova <- procD.lm(coords ~ species, 
                      data = Amb_gdf, iter = 999, 
                      RRPP = TRUE, print.progress = FALSE)
Amb_anova$aov.table

plot(Amb_anova, type = "diagnostics", outliers = TRUE)

# With size
Amb_anova_size <- procD.lm(coords ~ species*size, 
                           data = Amb_gdf, iter = 999, 
                           RRPP = TRUE, print.progress = FALSE)
Amb_anova_size$aov.table

plot(Amb_anova_size, type = "diagnostics", outliers = TRUE)

#Post-hoc comparisons

gp <-  interaction(Amb_gdf$species)
PW <- pairwise(Amb_anova, groups = gp, covariate = NULL)

summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)

## Phylogenetic signal ##

# Load in data #
require(phytools)

download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/GMM/Data/Amb_species?raw=true",   "Amb_species.txt")

# Read in tree

tree <- read.newick("Amb_species.txt")  #tree from Williams et al. 2013

par(mar=c(1,1,1,1))
tree$tip.label<-gsub("^", "A.", tree$tip.label)
plot(tree)

#Subset tree to include only GMM species
Amb_species<-unique(GMM_data_sub$species)
tips<-tree$tip.label
ii<-sapply(Amb_species,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
plotTree(tree,ftype="i")
#Tree did not include A.mavortium so I lumped that species with A.tigrinum
Amb_gdf$species<-gsub("A.mavortium", "A.tigrinum", Amb_gdf$species, fixed = TRUE)
Amb_gdf$species<-as.factor(Amb_gdf$species)

#Preformed a group PCA
library(Morpho)
gpca <- groupPCA(Amb_gdf$coords, Amb_gdf$species, rounds=0)
plot(gpca$groupmeans)

#Performed a Phylogenetic PCA based on group means
phylo.PCA <- gm.prcomp(gpca$groupmeans, phy = tree, align.to.phy = FALSE)
summary(phylo.PCA)

A_species<-attributes(gpca$groupmeans) #access attributes names
A_species<-(A_species$dimnames[[3]])
A_species<-as.factor(A_species)

#Plot phylogenetic PCA
plot(phylo.PCA, phylo = TRUE, main = "phylo PCA", col=A_species)

#3D plot of plylogenetic PCA
plotdat<-phylo.PCA$x[,1:3]
colnames(plotdat)<-c("","","")#prevent axis labels
obj<-phytools::phylomorphospace3d(tree,plotdat, method="dynamic",
                                  control=list(ftype="off",spin=FALSE, box=FALSE), cex.symbol=0.5)
spheres3d(phylo.PCA$x[,1:3], color=palette()[A_species], r=0.01)
bbox3d(color = c("white"),shininess=15, alpha=0.3,xat=c(10), xlab="x",yat=c(10), ylab="y",zat=c(10), zlab="z")
text3d((phylo.PCA$x[,1:3]+0.005), texts = substr(A_species,1,6))

#Test for phylogenetic signal, uses Blomberg’s K to test for strength and significance of phylogenetic signal.
physignal(gpca$groupmeans, tree, print.progress = F)

#Phylogenetic generalized least squares

avg_gdf<-geomorph.data.frame(coords=gpca$groupmeans, species=A_species) #make new geomorph dataframe with group mean coords

pgls<-procD.pgls(coords~species, phy=tree, data=avg_gdf, print.progress = F) #Phylogenetic generalized least squares
pgls$aov.table

#Compare evolutionary rates in different portions of the tree based on brownian motion

rate.comp<-compare.evol.rates(avg_gdf$coords, tree, gp=avg_gdf$species, method = c("permutation"), iter = 999, print.progress = F)

plot(rate.comp)
rate.comp$sigma.d.gp
rate.comp$pairwise.pvalue

## Discriminant Function Analysis ##

library(Morpho)

DFA<-CVA(GPA_landmarks_sub$coords, GMM_data_sub$species, cv = TRUE)

barplot(DFA$Var[,2]) # Variance explained by the canonical roots

# Assess the accuracy of jacknife #

accJack <- table(DFA$groups, DFA$class)
accJack
diag(prop.table(accJack, 1))
sum(accJack[row(accJack) == col(accJack)]) / sum(accJack)

# Plot first two DF axes #

DFA_cva <- data.frame(DFA$CVscores, species = DFA$groups)

ggplot(DFA_cva, aes(CV.1, CV.2)) +
  geom_point(aes(color = species)) + theme_classic()

#alternative plot
plot(DFA$CVscores, col=GMM_data_sub$species, pch=as.numeric(GMM_data_sub$species), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(DFA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(DFA$Var[2,2],1),"%")))
text(DFA$CVscores, as.character(GMM_data_sub$species), col=as.numeric(GMM_data_sub$species), cex=.7)

# Plot Mahalahobis distances as dendrogram #

dendroS=hclust(DFA$Dist$GroupdistMaha)
dendroS$labels=levels(GMM_data_sub$species)
par(mar=c(6.5,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="",
     ylab='Mahalahobis distance')


### K Nearest neighbor ###:Non-parametric
library(caret)
Atlas_PC_scores <- data.frame(Amb_PCA$x,species=GMM_data_sub$species)

set.seed(123)

KNNmodel <- train(
  species ~., data = Atlas_PC_scores, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center","scale"), #scale the data
  tuneLength = 20)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(Atlas_PC_scores[,1:17]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == Atlas_PC_scores$species) #overall accuracy

accKNN <- table(Atlas_PC_scores$species,predicted.classes)
accKNN
diag(prop.table(accKNN, 1))




