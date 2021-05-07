### Geometric Morphometrics of atlas ###

# Read in landmark tps file
library(curl)
library(geomorph)
f <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/GMM/Data/tps_compiled_5_4_21.TPS")
raw_data <- readland.tps(f, specID = c("imageID")) # the function "readland.tps" reads the landmark data in tps format and returns a 3D array of the coordinate data
plot(raw_data)
head(raw_data)

# Read in csv file of specimens
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/GMM/Data/Ambystoma_atlas_a_GMM_5_4_21.CSV")
specimenList <- read.csv(f2, header = FALSE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(specimenList)
library(dplyr)
species<-gsub("_.*","",specimenList$V2) #make a separate species vector
s <-cbind(specimenList, species)





## SUBSET DATA FOR FOSSILS ONLY ##

library(tidyverse)
new <- list(land = raw_data , species=s$species, specimen = s$V2)

# I first subset the original landmark data. 
# Data Processing

GPA_sub <- new 

# replace the "land" array (16 x 2 x 361 dimensions) in rawData with a list of 361 16 x 2 arrays
splitLand <- list(GPA_sub)
for (i in 1:155){
  
  splitLand[[i]] <-raw_data[1:8,1:2,i]
  
}
GPA_sub[["land"]] <- splitLand

# Specify a vector of the specimen names you want to extract
b <- data.frame(specimen = new$specimen, species = new$species) # creates matrix combining specimen numbers, centroid sizes, and species

aa <- dplyr::filter(b, grepl('41229', species))
aa <- aa$specimen
aa <- as.vector(aa)

FossilsToExtract <- c(aa)

# Create a tibble of just specimen names from your dataset and add a row ID column

Fossilspecimens <- tibble(specimenName=as.character(new$specimen))

Fossilspecimens <- rowid_to_column(Fossilspecimens, "ID")

# Filter that tibble to include only the ones you want to extract

Fossilspecimens <- filter(Fossilspecimens,specimenName %in% FossilsToExtract)

# Set up vectors of variables to pull out for each specimen

land <- vector("list",nrow(Fossilspecimens))

species <- vector()

specimen <- vector()

# Extract from rawData just the specimens for which you want data

j <- 0
for (i in Fossilspecimens$ID){
  
  j <- j + 1
  
  land[[j]] <- GPA_sub[["land"]][[i]]
  
  species <- c(species,as.character(GPA_sub[["species"]][[i]]))
  
  specimen <- c(specimen,as.character(GPA_sub[["specimen"]][[i]]))
  
}

# Assemble the extracted data into a format like that of your original dataset

GMM_data_fossil <- list("land"=unlist(land),"species"=as.factor(species),"specimen"=as.factor(specimen)) #raw data of fossil specimens


# the lines below recast the list in "land" to the original array format and attributes

dim(GMM_data_fossil$land) <- c(8,2,nrow(Fossilspecimens))

attributes(GMM_data_fossil$land)$dimnames[[3]] <- Fossilspecimens$specimenName





## SUBSET DATA TO REMOVE FOSSILS ##

# I first subset the original landmark data. 
# Data Processing

# Specify a vector of the specimen names you want to extract

noFossil<- dplyr::filter(b, !grepl('41229', species))
noFossil <- noFossil$specimen
noFossil <- as.vector(noFossil)

noFossilToExtract <- c(noFossil)

# Create a tibble of just specimen names from your dataset and add a row ID column

specimens <- tibble(specimenName=as.character(new$specimen))

specimens <- rowid_to_column(specimens, "ID")

# Filter that tibble to include only the ones you want to extract

specimens <- filter(specimens,specimenName %in% noFossilToExtract)

# Set up vectors of variables to pull out for each specimen

land <- vector("list",nrow(specimens))

species <- vector()

specimen <- vector()

# Extract from rawData just the specimens for which you want data

j <- 0
for (i in specimens$ID){
  
  j <- j + 1
  
  land[[j]] <- GPA_sub[["land"]][[i]]
  
  species <- c(species,as.character(GPA_sub[["species"]][[i]]))
  
  specimen <- c(specimen,as.character(GPA_sub[["specimen"]][[i]]))
  
}

# Assemble the extracted data into a format like that of your original dataset

GMM_data_noFossil <- list("land"=unlist(land),"species"=as.factor(species),"specimen"=as.factor(specimen)) #raw data of the smallest three juveniles of each species and the largest three adults of each species


# the lines below recast the list in "land" to the original array format and attributes

dim(GMM_data_noFossil$land) <- c(8,2,nrow(specimens))

attributes(GMM_data_noFossil$land)$dimnames[[3]] <- specimens$specimenName





## SUBSET DATA TO REMOVE A. subsalsum & A. ordinarium ##

# I first subset the original landmark data. 
# Data Processing

# Specify a vector of the specimen names you want to extract

Ambsub<- dplyr::filter(b, !grepl('41229|A.subsalsum|A.ordinarium', species))
Ambsub <- Ambsub$specimen
Ambsub <- as.vector(Ambsub)

AmbsubToExtract <- c(Ambsub)

# Create a tibble of just specimen names from your dataset and add a row ID column

specimens <- tibble(specimenName=as.character(new$specimen))

specimens <- rowid_to_column(specimens, "ID")

# Filter that tibble to include only the ones you want to extract

specimens <- filter(specimens,specimenName %in% AmbsubToExtract)

# Set up vectors of variables to pull out for each specimen

land <- vector("list",nrow(specimens))

species <- vector()

specimen <- vector()

# Extract from rawData just the specimens for which you want data

j <- 0
for (i in specimens$ID){
  
  j <- j + 1
  
  land[[j]] <- GPA_sub[["land"]][[i]]
  
  species <- c(species,as.character(GPA_sub[["species"]][[i]]))
  
  specimen <- c(specimen,as.character(GPA_sub[["specimen"]][[i]]))
  
}

# Assemble the extracted data into a format like that of your original dataset

GMM_data_sub <- list("land"=unlist(land),"species"=as.factor(species),"specimen"=as.factor(specimen)) #raw data of the smallest three juveniles of each species and the largest three adults of each species


# the lines below recast the list in "land" to the original array format and attributes

dim(GMM_data_sub$land) <- c(8,2,nrow(specimens))

attributes(GMM_data_sub$land)$dimnames[[3]] <- specimens$specimenName





# GENERALIZED PROCUSTES ANALYSIS ANALYSIS #

GPA_noFossil_landmarks <- gpagen(GMM_data_noFossil$land) # performs Generalized Procrustes analysis of landmarks and creates aligned Procrustes coordinates
GPA_sub_landmarks <- gpagen(GMM_data_sub$land)
GPA_Fossil_landmarks <- gpagen(GMM_data_fossil$land)

# CREATE DATAFRAMES 1) WITHOUT FOSSILS 2) WITHOUT FOSSILS & W/O A. SUBSALSUM AND A.ORDINARIUM 3) FOSSIL ONLY #
GMM_data_noFossil <-geomorph.data.frame(coords=GPA_noFossil_landmarks$coords,
                               size=GPA_noFossil_landmarks$Csize, species=GMM_data_noFossil$species)
GMM_data_noFossil$species <- factor(GMM_data_noFossil$species, levels = 
                                      c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                        "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder levels of species
GMM_data_sub <-geomorph.data.frame(coords=GPA_sub_landmarks$coords,
                                        size=GPA_sub_landmarks$Csize, species=GMM_data_sub$species)
GMM_data_sub$species <- factor(GMM_data_sub$species, levels = 
                                      c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                        "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.velasci")) # Reorder levels of species
GMM_data_fossil <-geomorph.data.frame(coords=GPA_Fossil_landmarks$coords,
                                        size=GPA_Fossil_landmarks$Csize, species=GMM_data_fossil$species)





## PRINCIPAL COMPONENT ANALYSIS ##

GMM_data_noFossil$coords <- two.d.array(GMM_data_noFossil$coords) #get the data in XY format for PCA
Amb_fossil_coords <- two.d.array(GMM_data_fossil$coords) #get the data in XY format for PCA

Amb_PCA <- prcomp(GMM_data_noFossil$coords) #PC analysis

# PROJECT FOSSIL DATA #
Amb_fossil_PCA <- predict(Amb_PCA, Amb_fossil_coords)
Fossil_PC_scores <- as.data.frame(Amb_fossil_PCA)
Fossil_PC_scores <- cbind(Fossil_PC_scores, species= GMM_data_fossil$species)

PC_scores <- as.data.frame(Amb_PCA$x)
PC_scores <- cbind(PC_scores, species= GMM_data_noFossil$species)
percentage <- round(Amb_PCA$sdev / sum(Amb_PCA$sdev) * 100, 2) # find percentage variance explained by PC's
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

All_PC_scores <- rbind(PC_scores, Fossil_PC_scores) # create a new dataframe with the original PC scores and the PC scores the fossils

speciescolors <- c("#666600", "#999900" ,"#99FFFF" ,"#FF9933" ,"#B266FF" ,"#1139BC", "#003BFD", "#00FD33", "#26B543", "#E50CF5", "#FF0909","#D71D1D", "#EC6A6A", "#8B0B0B", "#A54E4E")
#plot
library(ggplot2)
library(ggforce)
p<-ggplot(All_PC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species)) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_manual(name = "Species", values=c(speciescolors, "black", "black", "black", "black", "black")) + theme_classic()
p





### Statistical tests ###
# *removed A. subsalsum and A. ordinarium* see code above for removal process

## ANOVA ##

#sample sizes
Amb_sub<-data.frame(species=GMM_data_sub$species)
Amb_sub %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

# ANOVA Without Csize
Amb_anova <- procD.lm(coords ~ species, 
                      data = GMM_data_sub, iter = 1000, 
                      RRPP = TRUE, print.progress = FALSE)
Amb_anova$aov.table
plot(Amb_anova, type = "diagnostics", outliers = TRUE)

# ANOVA With Csize
Amb_anova_size <- procD.lm(coords ~ species*log(size), 
                           data = GMM_data_sub, iter = 1000, 
                           RRPP = TRUE, print.progress = FALSE)
Amb_anova_size$aov.table

plot(Amb_anova_size, type = "diagnostics", outliers = TRUE)

#Post-hoc pairwise comparisons on ANOVA w/o Csize

gp <-  interaction(GMM_data_sub$species)
PW <- pairwise(Amb_anova, groups = gp, covariate = NULL)
summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)




### K NEAREST NEIGHBOR   ###:Non-parametric

GMM_data_sub$coords <- two.d.array(GMM_data_sub$coords) #get the data in XY format for PCA
Amb_PCA_sub <- prcomp(GMM_data_sub$coords) #PCA

Amb_fossil_PCA2 <- predict(Amb_PCA_sub, Amb_fossil_coords) #PREDICT FOSSIL PC SCORES
Fossil_PC_scores2 <- as.data.frame(Amb_fossil_PCA2)

Atlas_PC_scores <- data.frame(Amb_PCA_sub$x,species=GMM_data_sub$species)

Atlas_PC_scores$ID <- seq.int(nrow(Atlas_PC_scores))

library(caret)
set.seed(123)

KNNmodel <- train(
  species ~., data = Atlas_PC_scores, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  tuneLength = 10)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(Atlas_PC_scores[,1:16]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == Atlas_PC_scores$species) #overall accuracy

accKNN <- table(Atlas_PC_scores$species,predicted.classes)
accKNN

# FOSSIL CLASSIFICATION #

library(class)
KnnTestPrediction_k13 <- knn(Atlas_PC_scores[,1:16], Fossil_PC_scores2,
                            Atlas_PC_scores$species, k=13, prob=TRUE)
KnnTestPrediction_k13






### RANDOM FOREST CLASSIFICATION ###:Non-parametric
library(randomForest)
set.seed(123)
Atlas.rf <- randomForest(species ~ ., data=Atlas_PC_scores, importance=TRUE)
print(Atlas.rf)
rf_acc <- Atlas.rf$confusion
rf_acc <- 1-rf_acc[,14] # percent correct classification
rf_acc

mean(Atlas.rf$predicted == Atlas_PC_scores$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred = predict(Atlas.rf, newdata = Fossil_PC_scores2[,1:16])
y_pred






# AMBYSTOMA CLADE CLASSIFICATION #

species=GMM_data_sub$species
clades <- recode(species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
Atlas_PC_scores_clade <- data.frame(Amb_PCA_sub$x,clades=clades)

#KNN#
set.seed(123)
KNNmodel_clades <- train(
  clades ~., data = Atlas_PC_scores_clade, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  tuneLength = 10)

plot(KNNmodel_clades) # plot accuracy vs k
KNNmodel_clades$bestTune # optimal k

predicted.clades <- KNNmodel_clades %>% predict(Atlas_PC_scores_clade[,1:16]) # predict class based on KNN model
head(predicted.clades)
mean(predicted.clades == Atlas_PC_scores_clade$clades) #overall accuracy

accKNN_clades <- table(Atlas_PC_scores_clade$clades,predicted.clades)
accKNN_clades

# FOSSIL CLASSIFICATION #

KnnTestPrediction_k9 <- knn(Atlas_PC_scores_clade[,1:16], Fossil_PC_scores2,
                            Atlas_PC_scores_clade$clades, k=9, prob=TRUE)
KnnTestPrediction_k9






# RANDOM FOREST WITH CLADES #
set.seed(123)
Atlas.rf_clades <- randomForest(clades ~ ., data=Atlas_PC_scores_clade, importance=TRUE)
print(Atlas.rf_clades)
rf_acc_clades <- Atlas.rf_clades$confusion
rf_acc_clades <- 1-rf_acc_clades[,9] # percent correct classification
rf_acc_clades

mean(Atlas.rf_clades$predicted == Atlas_PC_scores_clade$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_clade = predict(Atlas.rf_clades, newdata = Fossil_PC_scores2[,1:16])
y_pred_clade






# TEST WITH NUMBER OF SPECIMENS PER SPECIES EQUAL #

library(groupdata2)
library(dplyr)
library(ggplot2)
library(knitr) # kable()

set.seed(123) # For reproducibility

train_set <- balance(Atlas_PC_scores, 5, cat_col = "species") # 5 specimens per species
test_set <- Atlas_PC_scores[-train_set$ID, ] # test sample
train_set <- subset(train_set, select=-c(ID)) #remove ID column
test_set <- subset(test_set, select=-c(ID))


KnnTestPrediction_k5 <- knn(train_set[,1:16], test_set[,1:16],
                            train_set$species, k=5, prob=TRUE)

accKNN_euqal <- table(test_set$species,KnnTestPrediction_k5)
accKNN_euqal

mean(KnnTestPrediction_k5 == test_set$species) #overall accuracy







### MAHALAHOBIS DISTANCES AND NEIGHBOR JOINING ###

library(HDMD)
Mahala1 = pairwise.mahalanobis(All_PC_scores[,1:12], All_PC_scores$species, digits = 3)
names = rownames(Mahala1$means) #capture labels

mahala = sqrt(Mahala1$distance) #mahalanobis distance
rownames(mahala) = names #set rownames in the dissimilarity matrix
colnames(mahala) = names #set colnames in the dissimilarity matrix

mahala <- as.dist(mahala) #this is the mahalanobis dissimilarity matrix 

# dendroS <- hclust(mahala)

library("ape")
library(phytools)

# # plot(as.phylo(dendroS), type = "cladogram", cex = 0.6,
#      no.margin = TRUE)

tr <- nj(mahala) #neighbor joining

plot(unroot(as.phylo(tr)),type="unrooted",cex=0.6,
     use.edge.length=FALSE,lab4ut="axial",
     no.margin=TRUE)

plot((as.phylo(tr)),type="unrooted",cex=0.6,
     use.edge.length=TRUE,lab4ut="axial",
     no.margin=TRUE)



