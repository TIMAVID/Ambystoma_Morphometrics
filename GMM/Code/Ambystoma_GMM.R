## Load in data ##

download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/GMM/Data/GMM_data_noFossils.RData?raw=true",   "GMM_data_noFossil.RData")
load("GMM_data_noFossil.RData")

library(geomorph)


# Generalized procrustes analysis

GPA_landmarks <- gpagen(GMM_data_noFossil$land)

#Create geomorph data frame

Amb_gdf<-geomorph.data.frame(coords=GPA_landmarks$coords,
                             size=GPA_landmarks$Csize, species=GMM_data_noFossil$species)

## PCA ##

GPA_landmarks$coords <- two.d.array(GPA_landmarks$coords) #get the data in XY format for PCA
Amb_PCA <- prcomp(GPA_landmarks$coords)

### PCA vizualization ###

PC_scores <- as.data.frame(Amb_PCA$x)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(ggalt)
library(ggforce)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
percentage <- round(Amb_PCA$sdev / sum(Amb_PCA$sdev) * 100, 2)
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

# library(RColorBrewer)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# my.colors=rainbow(28) #set up color palette of rainbow colors with n = 14
# plot(1:28, pch=19, cex=2, col=my.colors)

GMM_data_noFossil$species <- factor(GMM_data_noFossil$species, levels = 
                                c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                  "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")) # Reorder species
levels(GMM_data_noFossil$species)

set.seed(123)
library(randomcoloR)
n <- 14
palette <- distinctColorPalette(n)


p<-ggplot(PC_scores,aes(x=PC1,y=PC2,color=GMM_data_noFossil$species ))
p<-p+geom_point(size =5)+theme + xlab(percentage[1]) + ylab(percentage[2]) +
  geom_encircle(expand=0, size = 3)+ theme_bw() + scale_color_manual(name = "Species", values=c("#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                                                              "#CDDADD", "#DC5956", "#E363BB")) + theme_classic()
# p + stat_ellipse()
p

### Load in Fossil data ###

download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/GMM/Data/GMM_data_fossil.RData?raw=true",   "GMM_data_fossil.RData")
load("GMM_data_fossil.RData")

# Generalized procrustes analysis

GPA_fossil_landmarks <- gpagen(GMM_data_fossil$land)

#Create geomorph data frame

Amb_fossil_gdf<-geomorph.data.frame(coords=GPA_fossil_landmarks$coords,
                             size=GPA_fossil_landmarks$Csize, species=GMM_data_fossil$species)
Amb_fossil_coords <- two.d.array(Amb_fossil_gdf$coords) #get the data in XY format for PCA

# Project data #
Amb_fossil_PCA <- predict(Amb_PCA, Amb_fossil_coords)
Fossil_PC_scores <- as.data.frame(Amb_fossil_PCA)

PC_scores <- cbind(PC_scores, genus= GMM_data_noFossil$species)
Fossil_PC_scores <- cbind(Fossil_PC_scores, genus= GMM_data_fossil$species)

All_PC_scores <- rbind(PC_scores, Fossil_PC_scores) # create a new dataframe with the original PC scores and the PC scores of your fossil
tail(All_PC_scores)
pointsToLabel <- as.character(GMM_data_fossil$species)

All_PC_scores$genus <- factor(All_PC_scores$genus, levels = 
                                        c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                          "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", pointsToLabel)) # Reorder species

levels(All_PC_scores$genus)


species <- c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
             "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")



pcaplot <- ggplot(data = All_PC_scores, mapping = aes(x = PC1, y = PC2,
                                           col = genus, label = genus)) # creates the initial plot with datapoints color-coded and unique symbols by each genus
pcaplot <- pcaplot + geom_encircle(expand=0, size = 2, data = All_PC_scores[!All_PC_scores$genus %in% pointsToLabel,])+ theme_bw()
pcaplot <- pcaplot + 
  geom_text(aes(PC1, PC2, label = genus), nudge_y = .003,
            check_overlap = FALSE, data = All_PC_scores[All_PC_scores$genus %in% pointsToLabel,])+ geom_point(data = All_PC_scores[All_PC_scores$genus %in% pointsToLabel,])

pcaplot <- pcaplot + scale_color_manual(name = "Species", breaks = c(species),
                                values=c("black", "black", "black", "black", "black", "#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                         "#CDDADD", "#DC5956", "#E363BB")) + theme_classic()
pcaplot

### Load in subset data ###

download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/GMM/Data/GMM_data_sub.RData?raw=true",   "GMM_data_sub.RData")
load("GMM_data_sub.RData")

# Generalized procrustes analysis

GMM_GPA_sub_coords <- gpagen(GMM_data_sub$land)

#Create geomorph data frame

Amb_gdf_sub<-geomorph.data.frame(coords=GMM_GPA_sub_coords$coords,
                             size=GMM_GPA_sub_coords$Csize, species=GMM_data_sub$species)

### Statistical tests ###

## ANOVA ##

#Check sample sizes

Amb_sub<-data.frame(species=GMM_data_sub$species)

t <- Amb_sub %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

write.table(t, file = "GMMSampleSize.txt", sep = ",", quote = FALSE, row.names = F)

# Without size
Amb_anova <- procD.lm(coords ~ species, 
                      data = Amb_gdf_sub, iter = 999, 
                      RRPP = TRUE, print.progress = FALSE)
Amb_anova$aov.table
?procD.lm
plot(Amb_anova, type = "diagnostics", outliers = TRUE)

# With size
Amb_anova_size <- procD.lm(coords ~ species*size, 
                           data = Amb_gdf_sub, iter = 999, 
                           RRPP = TRUE, print.progress = FALSE)
Amb_anova_size$aov.table

plot(Amb_anova_size, type = "diagnostics", outliers = TRUE)

?permudist

#Post-hoc comparisons

gp <-  interaction(Amb_gdf_sub$species)
PW <- pairwise(Amb_anova, groups = gp, covariate = NULL)
?pairwise
summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE)
t<- summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)
t <- t$pairwise.tables$P

write.table(t, file = "GMMPW.txt", sep = ",", quote = FALSE, row.names = T)

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
Amb_species<-unique(Amb_gdf_sub$species)
tips<-tree$tip.label
ii<-sapply(Amb_species,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
plotTree(tree,ftype="i")
#Tree did not include A.mavortium so I lumped that species with A.tigrinum
Amb_gdf_sub$species<-gsub("A.mavortium", "A.tigrinum", Amb_gdf_sub$species, fixed = TRUE)
Amb_gdf_sub$species<-as.factor(Amb_gdf_sub$species)

#Preformed a group PCA
library(Morpho)
gpca <- groupPCA(Amb_gdf_sub$coords, Amb_gdf_sub$species, rounds=0)
plot(gpca$groupmeans)
?groupPCA
#Performed a Phylogenetic PCA based on group means
phylo.PCA <- gm.prcomp(gpca$groupmeans, phy = tree, align.to.phy = FALSE)
summary(phylo.PCA)
?gm.prcomp
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
physignal(gpca$groupmeans, tree, print.progress = F, iter = 999)

#Phylogenetic generalized least squares

avg_gdf<-geomorph.data.frame(coords=gpca$groupmeans, species=A_species) #make new geomorph dataframe with group mean coords

pgls<-procD.pgls(coords~species, phy=tree, data=avg_gdf, print.progress = F, iter = 999) #Phylogenetic generalized least squares
pgls$aov.table
?procD.pgls

#Compare evolutionary rates in different portions of the tree based on brownian motion

rate.comp<-compare.evol.rates(avg_gdf$coords, tree, gp=avg_gdf$species, method = c("permutation"), iter = 999, print.progress = F)

plot(rate.comp)
rate.comp$sigma.d.gp
rate.comp$pairwise.pvalue

### K Nearest neighbor ###:Non-parametric

GMM_GPA_sub_coords$coords <- two.d.array(GMM_GPA_sub_coords$coords) #get the data in XY format for PCA
GPA_fossil_landmarks$coords <- two.d.array(GPA_fossil_landmarks$coords) #get the data in XY format for PCA
Amb_PCA_sub <- prcomp(GMM_GPA_sub_coords$coords)

Amb_fossil_PCA2 <- predict(Amb_PCA_sub, Amb_fossil_coords)
Fossil_PC_scores2 <- as.data.frame(Amb_fossil_PCA2)

library(caret)

Atlas_PC_scores <- data.frame(Amb_PCA_sub$x,species=GMM_data_sub$species)

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

t <- diag(prop.table(accKNN, 1))
write.table(t, file = "KNNAcc.txt", sep = ",", quote = FALSE, row.names = T)

# Fossil predictions #

library(class)
KnnTestPrediction_k7 <- knn(Atlas_PC_scores[,1:16], Fossil_PC_scores2,
                            Atlas_PC_scores$species, k=7, prob=TRUE)
t <- cbind(as.character(GMM_data_fossil$species), as.character(KnnTestPrediction_k7[1:5]), as.character(attr(KnnTestPrediction_k7, 'prob')))
write.table(t, file = "KNNFossil.txt", sep = ",", quote = FALSE, row.names = T)

KnnTestPrediction_k5 <- knn(Atlas_PC_scores[,1:16], Fossil_PC_scores2,
                            Atlas_PC_scores$species, k=5, prob=TRUE)
KnnTestPrediction_k5

KnnTestPrediction_k3 <- knn(Atlas_PC_scores[,1:16], Fossil_PC_scores2,
                            Atlas_PC_scores$species, k=3, prob=TRUE)
KnnTestPrediction_k3


## Discriminant Function Analysis ##

library(Morpho)

DFA<-CVA(GMM_GPA_sub_coords$coords, GMM_data_sub$species, cv = TRUE, rounds = 0)

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

# Predict fossils #

fossil_CVA_scores <- predict(DFA,GPA_fossil_landmarks$coords)

fossil_class <- classify(DFA, cv = FALSE, newdata = GPA_fossil_landmarks$coords)
fossil_class$class
fossil_class$posterior

#alternative plot
plot(DFA$CVscores, col=GMM_data_sub$species, pch=as.numeric(GMM_data_sub$species), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(DFA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(DFA$Var[2,2],1),"%")))
text(DFA$CVscores, as.character(GMM_data_sub$species), col=as.numeric(GMM_data_sub$species), cex=.7)
text(fossil_CVA_scores, as.character(GMM_data_fossil$species), cex=.7)
points(DFA$CVscores, col=as.numeric(palette))


### Plot Mahalahobis distances as dendrogram ###

df <- data.frame(All_PC_scores[,1:12], species = All_PC_scores$genus)
df <- df %>%
  group_by(species)

df <- aggregate(df[, 1:12], list(df$species), mean)

meltExpensesByMonth <- melt(df, id.vars=1)

library("tidyr")
library("magrittr")

meltExpensesByMonth <- pivot_wider(meltExpensesByMonth, 
            id_cols = variable,
            names_from = Group.1,
            values_from = c(value))

head(meltExpensesByMonth)

library(HDMD)
Mahala1 = pairwise.mahalanobis(All_PC_scores[,1:12], All_PC_scores$genus, digits = 3)
names = rownames(Mahala1$means) #capture labels

mahala = sqrt(Mahala1$distance) #mahalanobis distance
rownames(mahala) = names #set rownames in the dissimilarity matrix
colnames(mahala) = names #set colnames in the dissimilarity matrix

mahala <- as.dist(mahala) #this is the mahalanobis dissimilarity matrix 

dendroS <- hclust(mahala)

library("ape")

plot(as.phylo(dendroS), type = "unrooted", cex = 0.9,
     no.margin = TRUE)


