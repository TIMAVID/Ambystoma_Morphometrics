## Load in data ##
# Read in csv file of data
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/Linear_data/Data/Amb_linear_data_final.csv")
Amb_linear_data <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(Amb_linear_data)

### ATLAS ###

# Tidy data #
library(dplyr)

Atlas <-  Amb_linear_data[c(3, 8:13, 54:55)] #select only relevant atlas measurements

# REMOVE FOSSILS AND TUBERCULUM INTERGLENOIDEUM MEASUREMENT #
Atlas_wofossil <- dplyr::filter(Atlas, !grepl('41229*', species)) # remove fossils
row.names(Atlas_wofossil) <- Atlas_wofossil$specimen_num
Atlas_wofossil$species <- factor(Atlas_wofossil$species, levels = 
                                         c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                           "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels

Atlas_wofossil_noTub <- subset(Atlas_wofossil, select=-c(tub_interglen_extension)) # no tuberculum interglenoideum extension measurement
Atlas_wofossil_noTub <- na.omit(Atlas_wofossil_noTub) # remove rows with N/A's
Atlas_wofossil_noTub <- subset(Atlas_wofossil_noTub, select=-c(specimen_num))


## SUBSET DATA FOR FOSSILS ONLY ##

Atlas_fossil <- dplyr::filter(Atlas, grepl('41229*', species)) # fossils
Atlas_fossil <- na.omit(Atlas_fossil) # remove rows with N/A's
row.names(Atlas_fossil) <- Atlas_fossil$species
Atlas_fossil_noTub <- subset(Atlas_fossil, select=-c(tub_interglen_extension, specimen_num))




## PRINCIPAL COMPONENT ANALYSIS ##

Atlas.pca <- prcomp(Atlas_wofossil_noTub[c(1:6)], center = TRUE, scale = FALSE) # PCA
PC_scores <- as.data.frame(Atlas.pca$x)
PC_scores <- cbind(PC_scores, species= Atlas_wofossil_noTub$species)

percentage <- round(Atlas.pca$sdev / sum(Atlas.pca$sdev) * 100, 2)# find percentage variance explained by PC's
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

# PROJECT FOSSIL DATA #
Amb_fossil_PCA <- predict(Atlas.pca, Atlas_fossil_noTub[,1:6])
Fossil_PC_scores <- as.data.frame(Amb_fossil_PCA)
Fossil_PC_scores <- cbind(Fossil_PC_scores, species= Atlas_fossil_noTub$species)

All_PC_scores <- rbind(PC_scores, Fossil_PC_scores) # create a new dataframe with the original PC scores and the PC scores of your fossil
tail(All_PC_scores)

# PLOT #
speciescolors <- c("#666600", "#999900" ,"#99FFFF" ,"#FF9933" ,"#B266FF" ,"#1139BC", "#003BFD", "#00FD33", "#26B543", "#E50CF5", "#FF0909","#D71D1D", "#EC6A6A", "#8B0B0B", "#A54E4E")

library(ggplot2)
library(ggforce)
p<-ggplot(All_PC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species)) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(name = "Species", values=c(speciescolors, "black", "black", "black", "black", "black")) + theme_classic()
p


# TUBERCULUM INTERGLENOIDEUM PLOT #

library(EnvStats)
Tub_dat <- Atlas_wofossil[c(1,9)]
Tub_dat <- na.omit(Tub_dat) # remove rows with N/A's

ventral_extension_p <- ggplot(data = Tub_dat, aes(x = species, y = (tub_interglen_extension)))
ventral_extension_p <- ventral_extension_p + geom_boxplot(na.rm = TRUE)
ventral_extension_p <- ventral_extension_p + theme(axis.text.x = element_text(angle = 90))
ventral_extension_p <- ventral_extension_p + ylab("ventral extension (mm)") + stat_n_text() +theme_classic()
ventral_extension_p




### STATISTICAL TESTS ###
# *removed A. subsalsum and A. ordinarium* see code above for removal process
Atlas_wofossil_noTub_sub <- dplyr::filter(Atlas_wofossil_noTub, !grepl('A.subsalsum|A.ordinarium', species))
Atlas_wofossil_noTub_sub$species <- factor(Atlas_wofossil_noTub_sub$species, levels = 
                                   c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                     "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.velasci")) # Reorder species levels

# Sample sizes #
Atlas_wofossil_noTub_sub %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

## PERMUTATION MANOVA ##
library(RVAideMemoire)
require(vegan)

#Permutational Multivariate Analysis of Variance Using Euclidean Distance Matrices
adonis(Atlas_wofossil_noTub_sub[,1:6]~species,data=Atlas_wofossil_noTub_sub,method="mahalanobis") 

#pairwise comparisons between group levels with corrections for multiple testing
set.seed(123)
pairwise.perm.manova(Atlas_wofossil_noTub_sub[,1:6],Atlas_wofossil_noTub_sub$species,nperm=50, progress = FALSE) #needs more permutation but takes a long time
# t <- AtlasPPM$p.value
# t <-round(t, digits = 3)
# write.table(t, file = "Atlas linear PW", sep = ",", quote = FALSE, row.names = T)


# TUBERCULUM INTERGLENOIDEUM Permutational ANOVA #
set.seed(123)
perm.anova(Tub_dat$tub_interglen_extension ~ Tub_dat$species, nperm=1000)
#pairwise comparisons between group levels with corrections for multiple testing
library(rcompanion)
pairwisePermutationTest(tub_interglen_extension ~ species, data = Tub_dat, method = "fdr")

pairwise.perm.t.test(Tub_dat$tub_interglen_extension,Tub_dat$species,nperm=999,progress = FALSE)
# t <- t$p.value
# t <- round(t, 3)
# write.table(t, file = "Tub PW", sep = ",", quote = FALSE, row.names = T)





### K NEAREST NEIGHBOR   ###:Non-parametric
library(caret)

set.seed(123)
KNNmodel <- train(
  species ~., data = Atlas_wofossil_noTub_sub, method = "knn",
  trControl = trainControl("LOOCV", number =100),
  tuneLength = 10)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(Atlas_wofossil_noTub_sub[,1:6]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == Atlas_wofossil_noTub_sub$species) #overall accuracy

accKNN <- table(Atlas_wofossil_noTub_sub$species,predicted.classes)
accKNN
# t <- diag(prop.table(accKNN, 1))
# t <-round(t, digits = 2)
#write.table(t, file = "Atlas linear KNNAC", sep = ",", quote = FALSE, row.names = T)

# FOSSIL CLASSIFICATION #

library(class)
KnnTestPrediction_k15 <- knn(Atlas_wofossil_noTub_sub[,1:6], Atlas_fossil_noTub[,1:6],
                             Atlas_wofossil_noTub_sub$species, k=15, prob=TRUE)
KnnTestPrediction_k15





### RANDOM FOREST CLASSIFICATION ###:Non-parametric
library(randomForest)
set.seed(123)
Atlas.rf <- randomForest(species ~ ., data=Atlas_wofossil_noTub_sub, importance=TRUE)
print(Atlas.rf)
rf_acc <- Atlas.rf$confusion
rf_acc <- 1-rf_acc[,14] # percent correct classification
rf_acc

mean(Atlas.rf$predicted == Atlas_wofossil_noTub_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred = predict(Atlas.rf, newdata = Atlas_fossil_noTub[,1:6])
y_pred






# AMBYSTOMA CLADE CLASSIFICATION #

species <- Atlas_wofossil_noTub_sub$species
clades <- recode(species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
Atlas_wofossil_noTub_sub_clade <- data.frame(Atlas_wofossil_noTub_sub[,1:6],clades=clades)

#KNN#
set.seed(123)
KNNmodel_clades <- train(
  clades ~., data = Atlas_wofossil_noTub_sub_clade, method = "knn",
  trControl = trainControl("LOOCV", number =100),
  tuneLength = 10)

plot(KNNmodel_clades) # plot accuracy vs k
KNNmodel_clades$bestTune # optimal k

predicted.clades <- KNNmodel_clades %>% predict(Atlas_wofossil_noTub_sub_clade[,1:6]) # predict class based on KNN model
head(predicted.clades)
mean(predicted.clades == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy

accKNN_clades <- table(Atlas_wofossil_noTub_sub_clade$clades,predicted.clades)
accKNN_clades

# FOSSIL CLADE CLASSIFICATION #

KnnTestPrediction_k11 <- knn(Atlas_wofossil_noTub_sub_clade[,1:6], Atlas_fossil_noTub[,1:6],
                             Atlas_wofossil_noTub_sub_clade$clades, k=11, prob=TRUE)
KnnTestPrediction_k11






# RANDOM FOREST WITH CLADES #
set.seed(123)
Atlas.rf_clades <- randomForest(clades ~ ., data=Atlas_wofossil_noTub_sub_clade, importance=TRUE)
print(Atlas.rf_clades)
rf_acc_clades <- Atlas.rf_clades$confusion
rf_acc_clades <- 1-rf_acc_clades[,9] # percent correct classification
rf_acc_clades

mean(Atlas.rf_clades$predicted == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy

# FOSSIL CLADE CLASSIFICATION #
y_pred_clade = predict(Atlas.rf_clades, newdata = Atlas_fossil_noTub[,1:6])
y_pred_clade





### MAHALAHOBIS DISTANCES AND NEIGHBOR JOINING ###
Atlas_noNA <- na.omit(Atlas) # remove rows with N/A's
Atlas_noNA <- subset(Atlas_noNA, select=-c(specimen_num))

library(HDMD)
Mahala1 = pairwise.mahalanobis(Atlas_noNA[,1:7], Atlas_noNA$species, digits = 3)
names = rownames(Mahala1$means) #capture labels

mahala = sqrt(Mahala1$distance) #mahalanobis distance
rownames(mahala) = names #set rownames in the dissimilarity matrix
colnames(mahala) = names #set colnames in the dissimilarity matrix

mahala <- as.dist(mahala) #this is the mahalanobis dissimilarity matrix 

dendroS <- hclust(mahala)

library("ape")
library(phytools)

# plot(as.phylo(dendroS), type = "unrooted", cex = 0.6,lab4ut="axial",
     # no.margin = TRUE)

tr <- nj(mahala) #neighbor joining

plot(unroot(as.phylo(tr)),type="unrooted",cex=0.6,
     use.edge.length=FALSE,lab4ut="axial",
     no.margin=TRUE)

plot((as.phylo(tr)),type="unrooted",cex=0.6,
     use.edge.length=TRUE,lab4ut="axial",
     no.margin=TRUE)
# add.scale.bar(cex = 0.7, font = 2, col = "black")
# library(phytools)
# rr.interactive<-reroot(tr,interactive=TRUE)
# plotTree(rr.interactive)





# MEASUREMENT RELATIVE IMPORTANCE BASED ON RF #

library(tidyverse)
library(skimr)
library(knitr)
library(party)
library(GGally)
ggpairs(Atlas_wofossil_noTub_sub_clade[,1:6])

rf3 <- cforest(
  clades ~ .,
  data = Atlas_wofossil_noTub_sub_clade,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)

create_crfplot <- function(rf, conditional = TRUE){
  
  imp <- rf %>%
    varimp(conditional = conditional) %>% 
    as_tibble() %>% 
    rownames_to_column("Feature") %>% 
    rename(Importance = value)
  
  p <- ggplot(imp, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "#53cfff", width = 0.65) +
    coord_flip() + 
    theme_light(base_size = 20) +
    theme(axis.title.x = element_text(size = 15, color = "black"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = 15, color = "black"),
          axis.text.y  = element_text(size = 15, color = "black")) 
  return(p)
}

create_crfplot(rf3, conditional = TRUE)

# W/O HEIGHT #
set.seed(123)
Atlas.rf_clades_wo5 <- randomForest(clades ~ X1 + X2 + X3 + X4 + X6, data=Atlas_wofossil_noTub_sub_clade, importance=FALSE)
print(Atlas.rf_clades_wo5)
rf_acc_clades_wo5 <- Atlas.rf_clades_wo5$confusion
rf_acc_clades_wo5 <- 1-rf_acc_clades_wo5[,9] # percent correct classification
rf_acc_clades_wo5

mean(Atlas.rf_clades_wo5$predicted == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy

# W/O COTYLE HEIGHT #
set.seed(123)
Atlas.rf_clades_wo3 <- randomForest(clades ~ X1 + X2 + X4 + X5 + X6, data=Atlas_wofossil_noTub_sub_clade, importance=FALSE)
print(Atlas.rf_clades_wo3)
rf_acc_clades_wo3 <- Atlas.rf_clades_wo3$confusion
rf_acc_clades_wo3 <- 1-rf_acc_clades_wo3[,9] # percent correct classification
rf_acc_clades_wo3

mean(Atlas.rf_clades_wo3$predicted == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy

# W/O CENTRUM WIDTH #
set.seed(123)
Atlas.rf_clades_wo6 <- randomForest(clades ~ X1 + X2 +X3 + X4 + X5, data=Atlas_wofossil_noTub_sub_clade, importance=FALSE)
print(Atlas.rf_clades_wo6)
rf_acc_clades_wo6 <- Atlas.rf_clades_wo6$confusion
rf_acc_clades_wo6 <- 1-rf_acc_clades_wo6[,9] # percent correct classification
rf_acc_clades_wo6

mean(Atlas.rf_clades_wo6$predicted == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy

# W/O POSTERIOR END OF CENTRUM #
set.seed(123)
Atlas.rf_clades_wocentrum <- randomForest(clades ~ X2 +X3 + X4, data=Atlas_wofossil_noTub_sub_clade, importance=FALSE)
print(Atlas.rf_clades_wocentrum)
rf_acc_clades_wocen <- Atlas.rf_clades_wocentrum$confusion
rf_acc_clades_wocen <- 1-rf_acc_clades_wocen[,9] # percent correct classification
rf_acc_clades_wocen

mean(Atlas.rf_clades_wocentrum$predicted == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy


