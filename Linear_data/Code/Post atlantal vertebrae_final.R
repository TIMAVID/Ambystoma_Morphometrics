## Load in data ##
# Read in csv file of data---------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/Linear_data/Data/Amb_linear_data_final.csv")
Amb_linear_data <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(Amb_linear_data)

### Post atlantal vertebrae ###

# Tidy data #
library(dplyr)

PoAtVerts <-  Amb_linear_data[c(6:7, 14:48, 54:55)] #select only relevant Post atlantal measurements

# REMOVE FOSSILS #
PoAtVerts_wofossil <- dplyr::filter(PoAtVerts, !grepl('41229*', species)) # remove fossils
row.names(PoAtVerts_wofossil) <- PoAtVerts_wofossil$specimen_num
# PoAtVerts_wofossil <- subset(PoAtVerts_wofossil, select=-c(specimen_num))
PoAtVerts_wofossil <- droplevels(PoAtVerts_wofossil)
PoAtVerts_wofossil$species <- factor(PoAtVerts_wofossil$species, levels = 
                                            c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                              "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels


library(tidyverse)

# SUBSET DATA FOR EACH VERTEBRA #---------------------------------
T1 <- dplyr::select(PoAtVerts_wofossil, contains("a", ignore.case = FALSE),  contains("species")) %>% drop_na()


T4 <- dplyr::select(PoAtVerts_wofossil, contains("b", ignore.case = FALSE), contains("species")) %>% drop_na()

T8 <- dplyr::select(PoAtVerts_wofossil, contains("c", ignore.case = FALSE), contains("species")) %>% drop_na()
T8 <- subset(T8, select=-c(specimen_num))

T8_extension <- dplyr::select(PoAtVerts_wofossil, contains("Cen"), contains("species")) %>% drop_na()

T12<- dplyr::select(PoAtVerts_wofossil, contains("d", ignore.case = FALSE), contains("species")) %>% drop_na()

Sc <- dplyr::select(PoAtVerts_wofossil, num_range("X", 14:20), contains("species")) %>% drop_na()


# PCA PLOT OF MODERN SPECIMENS #---------------------------------

library(ggplot2)
library(ggforce)
library(ggfortify)
speciescolors <- c("#666600", "#C9D42D" ,"#42CEBC" ,"#F2AB1F" ,"#864ED0" ,"#261ACE", "#086AFD", "#08FD6A", "#0C8B3F", "#E50CF5", "#FF5E00","#FF0000", "#FF6A6A", "#D5930F", "#9E1616")


# T1 & T4 "Anterior trunk" #
#"Anterior"(1,4) comparative verts

T1_comb <- T1 %>%                          # Applying row_number function
  dplyr::mutate(Vert = "a")
T1_comb <- T1_comb %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)

T4_comb <- T4 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T4_comb <- T4_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "b")

TrunkAnt <- rbind(T1_comb, T4_comb)

T1_4.pca <- prcomp(TrunkAnt[c(1:7)], center = TRUE, scale = FALSE) # PCA

percentage_T1_4 <- round(T1_4.pca$sdev^2 / sum(T1_4.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_T1_4 <- paste(colnames(T1_4.pca$x),"(", paste( as.character(percentage_T1_4), "%", ")", sep="") )
T1_4.scores <- data.frame(T1_4.pca$x, species = TrunkAnt$species)
T1_4loadings <- data.frame(Variables = rownames(T1_4.pca$rotation), T1_4.pca$rotation)# Extract loadings of the variables

T1_4_plot<-ggplot(T1_4.scores,aes(x=PC1,y=PC2,color=species)) + 
  # geom_segment(data = T1loadings, aes(x = 0, y = 0, xend = (PC1),
  #                                      yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),
  #              color = "black") +annotate("text", x = (T1loadings$PC1), y = (T1loadings$PC2),
  #                                       label = T1loadings$Variables) +
  geom_point(size =2)+ xlab(percentage_T1_4[1]) + ylab(percentage_T1_4[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkAnt$species), values=c(speciescolors)) + 
  theme_classic() + ggtitle("T1 & T4") + theme(legend.position = "none")
T1_4_plot



# T8 "Middle trunk vertebrae #
T8.pca <- prcomp(T8[c(1:7)], center = TRUE, scale = FALSE) # PCA

percentage_T8 <- round(T8.pca$sdev^2 / sum(T8.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_T8 <- paste(colnames(T8.pca$x),"(", paste( as.character(percentage_T8), "%", ")", sep="") )
T8.scores <- data.frame(T8.pca$x, species = T8$species)
T8loadings <- data.frame(Variables = rownames(T8.pca$rotation), T8.pca$rotation)# Extract loadings of the variables

T8_plot<-ggplot(T8.scores,aes(x=PC1,y=PC2,color=species)) + 
  # geom_segment(data = T1loadings, aes(x = 0, y = 0, xend = (PC1),
  #                                      yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),
  #              color = "black") +annotate("text", x = (T1loadings$PC1), y = (T1loadings$PC2),
  #                                       label = T1loadings$Variables) +
  geom_point(size =2)+ xlab(percentage_T8[1]) + ylab(percentage_T8[2]) +
  scale_color_manual(name = "Species", breaks=levels(T8$species), values=c(speciescolors)) + 
  theme_classic() + ggtitle("T8") + theme(legend.position = "none")
T8_plot


# T12 "Posterior trunk vertebrae #
T12.pca <- prcomp(T12[c(1:7)], center = TRUE, scale = FALSE) # PCA

percentage_T12 <- round(T12.pca$sdev^2 / sum(T12.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_T12 <- paste(colnames(T12.pca$x),"(", paste(as.character(percentage_T12), "%", ")", sep="") )
T12.scores <- data.frame(T12.pca$x, species = T8$species)
T12loadings <- data.frame(Variables = rownames(T12.pca$rotation), T12.pca$rotation)# Extract loadings of the variables

T12_plot<-ggplot(T12.scores,aes(x=PC1,y=PC2,color=species)) + 
  # geom_segment(data = T1loadings, aes(x = 0, y = 0, xend = (PC1),
  #                                      yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),
  #              color = "black") +annotate("text", x = (T1loadings$PC1), y = (T1loadings$PC2),
  #                                       label = T1loadings$Variables) +
  geom_point(size =2)+ xlab(percentage_T12[1]) + ylab(percentage_T12[2]) +
  scale_color_manual(name = "Species", breaks=levels(T12$species), values=c(speciescolors)) + 
  theme_classic() + ggtitle("T12") + theme(legend.position = "none")
T12_plot


# SC Sacral vertebrae #
Sc.pca <- prcomp(Sc[c(1:7)], center = TRUE, scale = FALSE) # PCA

percentage_SC <- round(Sc.pca$sdev^2 / sum(Sc.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_SC <- paste(colnames(Sc.pca$x),"(", paste( as.character(percentage_SC), "%", ")", sep="") )
SC.scores <- data.frame(Sc.pca$x, species = Sc$species)
SCloadings <- data.frame(Variables = rownames(Sc.pca$rotation), Sc.pca$rotation)# Extract loadings of the variables

SC_plot<-ggplot(SC.scores,aes(x=PC1,y=PC2,color=species)) + 
  # geom_segment(data = SCloadings, aes(x = 0, y = 0, xend = (PC1),
  #                                      yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")),
  #              color = "black") +annotate("text", x = (T1loadings$PC1), y = (T1loadings$PC2),
  #                                       label = T1loadings$Variables) +
  geom_point(size =2)+ xlab(percentage_SC[1]) + ylab(percentage_SC[2]) +
  scale_color_manual(name = "Species", breaks=levels(Sc$species), values=c(speciescolors)) + 
  theme_classic() + ggtitle("SC") + theme(legend.position = "none")
SC_plot


# PLOT ALL TOGETHER #
par(oma=c(0,0,2,0))
gridExtra::grid.arrange(T1_4_plot, T8_plot ,T12_plot,SC_plot,nrow = 2)




# *removed A. subsalsum and A. ordinarium* #---------------------------------

#Anterior vertebrae

TrunkAnt_sub <- dplyr::filter(TrunkAnt, !grepl('A.subsalsum|A.ordinarium', species))
TrunkAnt_sub$species <- droplevels(TrunkAnt_sub$species)

#Mid vertebrae
Trunk8_sub <- dplyr::filter(T8, !grepl('A.subsalsum|A.ordinarium', species))
Trunk8_sub$species <- droplevels(Trunk8_sub$species)


#Posterior vertebrae
Trunk12_sub <- dplyr::filter(T12, !grepl('A.subsalsum|A.ordinarium', species))
Trunk12_sub$species <- droplevels(Trunk12_sub$species)


# Sacral vertebrae
Sc_sub <- dplyr::filter(Sc, !grepl('A.subsalsum|A.ordinarium', species))
Sc_sub$species <- droplevels(Sc_sub$species)







### K NEAREST NEIGHBOR   ###:Non-parametric---------------------------------
library(caret)

#make KNN model using LOOCV to find optimal k

# T1&4 VERTEBRAE #

# LOOCV WITH REPLICATION
library(foreach)
library(doParallel)
ncore <- detectCores()
registerDoParallel(cores=ncore)

set.seed(123)
runs <- 100
system.time({
  fishT1 <- foreach(icount(runs)) %dopar% {
    caret::train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = caret::trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = TrunkAnt_sub)$results
  }
}) #repeated KNN model using LOOCV to find optimal k


fishT1 <- map_dfr(fishT1,`[`, c("k", "Accuracy", "Kappa"))
kfishT1 <- fishT1 %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishT1 <- kfishT1[1,] # k with highest accuracy


set.seed(123)
predicted.classes_T1 <- train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                              method     = "knn",
                              tuneGrid   = expand.grid(k = 3),
                              trControl  = trainControl(method  = "LOOCV"),
                              metric     = "Accuracy",
                              data       = TrunkAnt_sub)$pred # predict class based on KNN model
mean(predicted.classes_T1$pred == predicted.classes_T1$obs) #overall accuracy

accKNNT1 <- table(predicted.classes_T1$obs,predicted.classes_T1$pred)
accKNNT1
t <- diag(prop.table(accKNNT1, 1))
t <-round(t, digits = 2)
write.table(t, file = "T1 species KNNAC", sep = ",", quote = FALSE, row.names = T)





## T8 VERTEBRAE ##

# LOOCV WITH REPLICATION


set.seed(123)
runs <- 100
system.time({
  fishT8 <- foreach(icount(runs)) %dopar% {
    train(species ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = Trunk8_sub)$results
  }
}) #repeated KNN model using LOOCV to find optimal k

fishT8 <- map_dfr(fishT8,`[`, c("k", "Accuracy", "Kappa"))
kfishT8 <- fishT8 %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishT8 <- kfishT8[1,] # k with highest accuracy

set.seed(123)
predicted.classes_T8 <- train(species ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c,
                              method     = "knn",
                              tuneGrid   = expand.grid(k = 4),
                              trControl  = trainControl(method  = "LOOCV"),
                              metric     = "Accuracy",
                              data       = Trunk8_sub)$pred # predict class based on KNN model
mean(predicted.classes_T8$pred == predicted.classes_T8$obs) #overall accuracy

accKNNT8 <- table(predicted.classes_T8$obs,predicted.classes_T8$pred)
accKNNT8
# t <- diag(prop.table(accKNNT8, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "T8 species KNNAC", sep = ",", quote = FALSE, row.names = T)



## T12 VERTEBRAE ##

# LOOCV WITH REPLICATION

set.seed(123)
runs <- 100
system.time({
  fishT12 <- foreach(icount(runs)) %dopar% {
    train(species ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = Trunk12_sub)$results
  }
}) #repeated KNN model using LOOCV to find optimal k

fishT12 <- map_dfr(fishT12,`[`, c("k", "Accuracy", "Kappa"))
kfishT12 <- fishT12 %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishT12 <- kfishT12[1,] # k with highest accuracy

set.seed(123)
predicted.classes_T12 <- train(species ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d,
                              method     = "knn",
                              tuneGrid   = expand.grid(k = 4),
                              trControl  = trainControl(method  = "LOOCV"),
                              metric     = "Accuracy",
                              data       = Trunk12_sub)$pred # predict class based on KNN model
mean(predicted.classes_T12$pred == predicted.classes_T12$obs) #overall accuracy

accKNNT12 <- table(predicted.classes_T12$obs,predicted.classes_T12$pred)
accKNNT12
# t <- diag(prop.table(accKNNT12, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "T12 species KNNAC", sep = ",", quote = FALSE, row.names = T)


### RANDOM FOREST CLASSIFICATION ###:Non-parametric---------------------------------
library(randomForest)

# ANTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_ant <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkAnt_sub, importance=FALSE)
print(Atlas.rf_ant)
rf_acc_ant <- Atlas.rf_ant$confusion
rf_acc_ant <- 1-rf_acc_ant[,14] # percent correct classification
rf_acc_ant

# t <- rf_acc_ant
# t <-round(t, digits = 2)
# write.table(t, file = "Anterior verts RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_ant$predicted == TrunkAnt_sub$species) #overall accuracy


# T8 #
set.seed(123)
Atlas.rf_T8 <- randomForest(species ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c, data=Trunk8_sub, importance=FALSE)
print(Atlas.rf_T8)
rf_acc_T8 <- Atlas.rf_T8$confusion
rf_acc_T8 <- 1-rf_acc_T8[,14] # percent correct classification
rf_acc_T8

# t <- rf_acc_T8
# t <-round(t, digits = 2)
# write.table(t, file = "T8 RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_T8$predicted == Trunk8_sub$species) #overall accuracy


# T12 #
set.seed(123)
Atlas.rf_T12 <- randomForest(species ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d, data=Trunk12_sub, importance=FALSE)
print(Atlas.rf_T12)
rf_acc_T12 <- Atlas.rf_T12$confusion
rf_acc_T12 <- 1-rf_acc_T12[,14] # percent correct classification
rf_acc_T12

# t <- rf_acc_T12
# t <-round(t, digits = 2)
# write.table(t, file = "T12 RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_T12$predicted == Trunk12_sub$species) #overall accuracy


# AMBYSTOMA CLADE CLASSIFICATION #---------------------------------
TrunkAnt_sub$clades <- dplyr::recode(TrunkAnt_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
Trunk8_sub$clades <- dplyr::recode(Trunk8_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
Trunk12_sub$clades <- dplyr::recode(Trunk12_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
Sc_sub$clades <- dplyr::recode(Sc_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")

### K NEAREST NEIGHBOR CLADES ###:Non-parametric---------------------------------
library(caret)

## T8 VERTEBRAE ##

# LOOCV WITH REPLICATION


set.seed(123)
runs <- 100
system.time({
  fishT8clades <- foreach(icount(runs)) %dopar% {
    train(clades ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = Trunk8_sub)$results
  }
}) #repeated KNN model using LOOCV to find optimal k

fishT8clades <- map_dfr(fishT8clades,`[`, c("k", "Accuracy", "Kappa"))
kfishT8clades <- fishT8clades %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishT8clades <- kfishT8clades[1,] # k with highest accuracy

set.seed(123)
predicted.classes_T8clades <- train(clades ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c,
                              method     = "knn",
                              tuneGrid   = expand.grid(k = 6),
                              trControl  = trainControl(method  = "LOOCV"),
                              metric     = "Accuracy",
                              data       = Trunk8_sub)$pred # predict class based on KNN model
mean(predicted.classes_T8clades$pred == predicted.classes_T8clades$obs) #overall accuracy

accKNNT8clades <- table(predicted.classes_T8clades$obs,predicted.classes_T8clades$pred)
accKNNT8clades
t <- diag(prop.table(accKNNT8clades, 1))
t <-round(t, digits = 2)
write.table(t, file = "T8 clades KNNAC", sep = ",", quote = FALSE, row.names = T)



## T12 VERTEBRAE ##

# LOOCV WITH REPLICATION

set.seed(123)
runs <- 100
system.time({
  fishT12clades <- foreach(icount(runs)) %dopar% {
    train(clades ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = Trunk12_sub)$results
  }
}) #repeated KNN model using LOOCV to find optimal k

fishT12clades <- map_dfr(fishT12clades,`[`, c("k", "Accuracy", "Kappa"))
kfishT12clades <- fishT12clades %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishT12clades <- kfishT12clades[1,] # k with highest accuracy

set.seed(123)
predicted.classes_T12clades <- train(clades ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d,
                               method     = "knn",
                               tuneGrid   = expand.grid(k = 7),
                               trControl  = trainControl(method  = "LOOCV"),
                               metric     = "Accuracy",
                               data       = Trunk12_sub)$pred # predict class based on KNN model
mean(predicted.classes_T12clades$pred == predicted.classes_T12clades$obs) #overall accuracy

accKNNT12clades <- table(predicted.classes_T12clades$obs,predicted.classes_T12clades$pred)
accKNNT12clades
t <- diag(prop.table(accKNNT12clades, 1))
t <-round(t, digits = 2)
write.table(t, file = "T12 clades KNNAC", sep = ",", quote = FALSE, row.names = T)



### RANDOM FOREST CLADE CLASSIFICATION ###:Non-parametric---------------------------------
library(randomForest)

# T8 #
set.seed(123)
Atlas.rf_T8clades <- randomForest(clades ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c, data=Trunk8_sub, importance=FALSE)
print(Atlas.rf_T8clades)
rf_acc_T8clades <- Atlas.rf_T8clades$confusion
rf_acc_T8clades <- 1-rf_acc_T8clades[,9] # percent correct classification
rf_acc_T8clades

t <- rf_acc_T8clades
t <-round(t, digits = 2)
write.table(t, file = "T8 RFAC clades", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_T8clades$predicted == Trunk8_sub$clades) #overall accuracy


# T12 #
set.seed(123)
Atlas.rf_T12clades <- randomForest(clades ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d, data=Trunk12_sub, importance=FALSE)
print(Atlas.rf_T12clades)
rf_acc_T12clades <- Atlas.rf_T12clades$confusion
rf_acc_T12clades <- 1-rf_acc_T12clades[,9] # percent correct classification
rf_acc_T12clades

t <- rf_acc_T12clades
t <-round(t, digits = 2)
write.table(t, file = "T12 RFAC clades", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_T12clades$predicted == Trunk12_sub$clades) #overall accuracy





# MEASUREMENT RELATIVE IMPORTANCE BASED ON RF #---------------------------------
library(tidyverse)
library(skimr)
library(knitr)
library(party)
library(GGally)




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
          axis.text.y  = element_text(size = 15, color = "black"))  + theme_classic()
  return(p)
}

# CONDITIONAL PERUMATATION IMPORTANCE # *long time to run---------------------------------

Atlas.rf_ant_imp <- cforest(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, 
  data=TrunkAnt_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
Atlas.rf_ant_impplot <- create_crfplot(Atlas.rf_ant_imp, conditional = TRUE)


Atlas.rf_mid_imp <- cforest(
  clades ~ X7c + X8c + X9c + X10c + X11c + X12c + X13c, data=Trunk8_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
Atlas.rf_mid_impplot <- create_crfplot(Atlas.rf_mid_imp, conditional = TRUE)


Atlas.rf_post_imp <- cforest(
  clades ~ X7d + X8d + X9d + X10d + X11d + X12d + X13d, data=Trunk12_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
Atlas.rf_post_impplot <- create_crfplot(Atlas.rf_post_imp, conditional = TRUE)


Atlas.rf_sc_imp <- cforest(
  clades ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data=Sc_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
Atlas.rf_sc_impplot <- create_crfplot(Atlas.rf_sc_imp, conditional = TRUE)


# PLOT ALL TOGETHER #
par(oma=c(0,0,2,0))
gridExtra::grid.arrange(Atlas.rf_ant_impplot, Atlas.rf_mid_impplot, Atlas.rf_post_impplot, Atlas.rf_sc_impplot,nrow = 2)


#With important VH variable removed

library(randomForest)
# T8 #
set.seed(123)
Atlas.rf_T8clades <- randomForest(clades ~ X7c + X8c + X9c + X10c + X11c + X12c, data=Trunk8_sub, importance=FALSE)
print(Atlas.rf_T8clades)
rf_acc_T8clades <- Atlas.rf_T8clades$confusion
rf_acc_T8clades <- 1-rf_acc_T8clades[,9] # percent correct classification
rf_acc_T8clades

mean(Atlas.rf_T8clades$predicted == Trunk8_sub$clades) #overall accuracy

# T12 #
set.seed(123)
Atlas.rf_T12clades <- randomForest(clades ~ X7d + X8d + X9d + X10d + X11d + X12d, data=Trunk12_sub, importance=FALSE)
print(Atlas.rf_T12clades)
rf_acc_T12clades <- Atlas.rf_T12clades$confusion
rf_acc_T12clades <- 1-rf_acc_T12clades[,9] # percent correct classification
rf_acc_T12clades

mean(Atlas.rf_T12clades$predicted == Trunk12_sub$clades) #overall accuracy


# Sacral vertebra #
set.seed(123)
Atlas.rf_sc_clade <- randomForest(clades ~ X14 + X15 + X16 + X17 + X18 + X19, data=Sc_sub, importance=FALSE)
print(Atlas.rf_sc_clade)
rf_acc_sc_clade <- Atlas.rf_sc_clade$confusion
rf_acc_sc_clade <- 1-rf_acc_sc_clade[,9] # percent correct classification
rf_acc_sc_clade

mean(Atlas.rf_sc_clade$predicted == Sc_sub$clades) #overall accuracy



### CALCULATE ZYGAPOPHYSEAL RATIOS ###---------------------------------
library(tibble)
zygapophyses <- dplyr::select(PoAtVerts_wofossil, contains("X10"), contains("X11"), contains("X12"),contains("X17"), contains("X18"), contains("X19"), contains("species"), contains("specimen"))


# T1
ratio1 <- zygapophyses %>%                          
  dplyr::select(contains("a",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Zyratioa = (X10a+X11a)/X12a) %>% 
  dplyr::mutate(Vert = "a") %>% drop_na()
ratio1$specimen_num <- as.character(ratio1$specimen_num)

#T4
ratio4 <- zygapophyses %>%                          
  dplyr::select(contains("b",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Zyratiob = (X10b+X11b)/X12b) %>% 
  dplyr::mutate(Vert = "b") %>% drop_na()
ratio4$specimen_num <- as.character(ratio4$specimen_num)

#T8
ratio8 <- zygapophyses %>%                        
  dplyr::select(contains("c",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Zyratioc = (X10c+X11c)/X12c) %>% 
  dplyr::mutate(Vert = "c") %>% drop_na()
ratio8$specimen_num <- as.character(ratio8$specimen_num)

#T12
ratio12 <- zygapophyses %>%                        
  dplyr::select(contains("d",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Zyratiod = (X10d+X11d)/X12d) %>% 
  dplyr::mutate(Vert = "d") %>% drop_na()
ratio12$specimen_num <- as.character(ratio12$specimen_num)

#SC
ratioSC <- zygapophyses %>%                        
  dplyr::select(contains("X17"), contains("X18"), contains("X19"), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(ZyratioSc = (X17+X18)/X19) %>% 
  dplyr::mutate(Vert = "Sc") %>% drop_na()
ratioSC$specimen_num <- as.character(ratioSC$specimen_num)


library(plyr)
ZYRATIO <- join_all(list(ratio1, ratio4,ratio8,ratio12,ratioSC), by='specimen_num', type='left')

library("reshape2")
ZYRATIO_data_long <- melt(ZYRATIO, id=c("specimen_num", "species", "Vert"))  # convert to long format
ZYRATIO_data_long <- na.omit(ZYRATIO_data_long)

ZySummary <- ZYRATIO_data_long %>%
  dplyr::group_by(species, variable) %>%
  dplyr::summarise(mean = mean(value), min = min(value), max = max(value)) %>% 
  dplyr::filter(grepl("Zy*", variable))
ZySummary$species <- factor(ZySummary$species, levels = 
                                       c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                         "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels

speciescolors <- c("#666600", "#C9D42D" ,"#42CEBC" ,"#F2AB1F" ,"#864ED0" ,"#261ACE", "#086AFD", "#08FD6A", "#0C8B3F", "#E50CF5", "#FF5E00","#FF0000", "#FF6A6A", "#D5930F", "#9E1616")

ZyPlot <- ggplot(ZySummary, aes(x=variable, y=mean, group=species, color=species)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.1, position=position_dodge(0.15)) +
  geom_line() + geom_point()+
  scale_color_manual(name = "Species", values=c(speciescolors)) + theme_minimal() + xlab("Vertebra") + ylab("Zygapoheseal ratio") 
ZyPlot <- ZyPlot + facet_wrap(~species, ncol = 15)  # wrap data 'by' family into 4 columns
ZyPlot + scale_x_discrete(breaks=c("Zyratioa","Zyratiob","Zyratioc", "Zyratiod", "ZyratioSc"),
                           labels=c("1", "4", "8", "12", "Sc")) + theme(legend.position = "bottom")




### CALCULATE CENTRUM RATIOS ###---------------------------------

centrum <- dplyr::select(PoAtVerts_wofossil, contains("X7"), contains("X8"), contains("X9"),contains("X14"), contains("X15"), contains("X16"), contains("species"), contains("specimen"))


# T1
centrum_ratio1 <- centrum %>%                          
  dplyr::select(contains("a",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Cenratioa = (X7a)/X8a) %>% 
  dplyr::mutate(Vert = "a") %>% drop_na()
centrum_ratio1$specimen_num <- as.character(centrum_ratio1$specimen_num)

#T4
centrum_ratio4 <- centrum %>%                          
  dplyr::select(contains("b",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Cenratiob = (X7b)/X8b) %>% 
  dplyr::mutate(Vert = "b") %>% drop_na()
centrum_ratio4$specimen_num <- as.character(centrum_ratio4$specimen_num)

#T8
centrum_ratio8 <- centrum %>%                        
  dplyr::select(contains("c",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Cenratioc = (X7c)/X8c) %>% 
  dplyr::mutate(Vert = "c") %>% drop_na()
centrum_ratio8$specimen_num <- as.character(centrum_ratio8$specimen_num)

#T12
centrum_ratio12 <- centrum %>%                        
  dplyr::select(contains("d",ignore.case = FALSE), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(Cenratiod = (X7d)/X8d) %>% 
  dplyr::mutate(Vert = "d") %>% drop_na()
centrum_ratio12$specimen_num <- as.character(centrum_ratio12$specimen_num)

#SC
centrum_ratioSC <- centrum %>%                        
  dplyr::select(contains("X14"), contains("X15"), contains("X16"), contains("species"), contains("specimen")) %>% 
  dplyr::mutate(CenratioSc = (X14)/X15) %>% 
  dplyr::mutate(Vert = "Sc") %>% drop_na()
centrum_ratioSC$specimen_num <- as.character(centrum_ratioSC$specimen_num)


library(plyr)
CENRATIO <- join_all(list(centrum_ratio1, centrum_ratio4,centrum_ratio8,centrum_ratio12,centrum_ratioSC), by='specimen_num', type='left')

library("reshape2")
CENRATIO_data_long <- melt(CENRATIO, id=c("specimen_num", "species", "Vert"))  # convert to long format
CENRATIO_data_long <- na.omit(CENRATIO_data_long)

CENSummary <- CENRATIO_data_long %>%
  dplyr::group_by(species, variable) %>%
  dplyr::summarise(mean = mean(value), min = min(value), max = max(value)) %>% 
  dplyr::filter(grepl("Cen*", variable))
CENSummary$species <- factor(CENSummary$species, levels = 
                              c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels

CenPlot <- ggplot(CENSummary, aes(x=variable, y=mean, group=species, color=species)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.1, position=position_dodge(0.15)) +
  geom_line() + geom_point()+
  scale_color_manual(name = "Species", values=c(speciescolors)) + theme_minimal() + xlab("Vertebra") + ylab("Centrum ratio") 
CenPlot <- CenPlot + facet_wrap(~species, ncol = 15)  # wrap data 'by' family into 4 columns
CenPlot + scale_x_discrete(breaks=c("Cenratioa","Cenratiob","Cenratioc", "Cenratiod", "CenratioSc"),
                           labels=c("1", "4", "8", "12", "Sc")) + theme(legend.position = "bottom")




# T8 POSTERIOR EXTENSION MEASUREMENTS PLOTS #---------------------------------

library(EnvStats)
T8_extension

T8_extension<-dplyr::mutate(T8_extension ,ratio = Cen_to_NeuAr/Cen_to_PoZy, .before = Cen_to_PoZy) #add new column for ratio

library(tidyr)
data_long <- gather(T8_extension, Type , Measurement, Cen_to_NeuAr:Cen_to_PoZy, factor_key=TRUE) #convert to long format

library(ggplot2)
s <- ggplot(data_long, aes(Type, Measurement, fill = Type)) + 
  geom_boxplot(position = "dodge") + theme_classic() + ylab("Value") +
  theme(legend.position="bottom") + scale_fill_discrete(name = "Measurement Type", labels= c("NSPE", "NSPE/PoZyPE", "PoZyPE"))
s <- s + facet_wrap(~species, ncol = 5) + stat_n_text()
s


# ## Phylogenetic signal ##---------------------------------
# 
# # Load in data #
# require(phytools)
# 
# f3 <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/GMM/Data/Amb_species")
# 
# # Read in tree
# 
# tree <- read.newick(f3)  #tree from Williams et al. 2013
# 
# par(mar=c(1,1,1,1))
# tree$tip.label<-gsub("^", "A.", tree$tip.label)
# plot(tree)
# 
# #Subset tree to include only GMM species
# Amb_species<-unique(T8$species)
# tips<-tree$tip.label
# ii<-sapply(Amb_species,function(x,y) grep(x,y)[1],y=tips)
# tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
# plotTree(tree,ftype="i")
# 
# #Tree did not include A.mavortium so I lumped that species with A.tigrinum
# T8_sub_tig<- ZySummary
# T8_sub_tig$species<-gsub("A.mavortium", "A.tigrinum", T8_sub_tig$species, fixed = TRUE)
# T8_sub_tig$species<-gsub("A.subsalsum", "A.ordinarium", T8_sub_tig$species, fixed = TRUE)
# T8_sub_tig$species<-gsub("A.velasci", "A.ordinarium", T8_sub_tig$species, fixed = TRUE)
# T8_sub_tig$species<-as.factor(T8_sub_tig$species)
# T8_sub_tig$species <- factor(T8_sub_tig$species, levels = 
#                                c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
#                                  "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.ordinarium")) # Reorder species
# T8_sub_tig <- T8_sub_tig %>% 
#   dplyr::filter(grepl("Zyratioc", variable)) %>% 
#   dplyr::group_by(species) %>%
#   dplyr::summarise(Zymean = mean(mean))
# 
# library(geomorph)
# test <- (T8_sub_tig$mean)
# names <- T8_sub_tig$species
# setNames(test, T8_sub_tig$species)
# 
# test$dimnames[[1]] <- as.character(names)
# 
# physignal(T8_sub_tig$mean, tree, print.progress = F, iter = 999)
# 
# 
# names(T8_sub_tig$mean)














# MAKE CONSERVATIVE REGIONAL VERTEBRAE GROUPING SUBDATASETS "ANTERIOR, MIDDLE, POSTERIOR" #---------------------------------

#"Anterior"(1,4) comparative verts

T1_comb <- T1 %>%                          # Applying row_number function
  dplyr::mutate(Vert = "a")
T1_comb <- T1_comb %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)

T4_comb <- T4 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T4_comb <- T4_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "b")

TrunkAnt <- rbind(T1_comb, T4_comb)
TrunkAnt.pca <- prcomp(TrunkAnt[c(1:7)], center = TRUE, scale = FALSE) # PCA

#"Mid"(4,8,12) comparative verts
T8_comb <- T8 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T8_comb <- T8_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "c")

T12_comb <- T12 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T12_comb <- T12_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "d")

TrunkMid <- rbind(T4_comb, T8_comb, T12_comb)
TrunkMid.pca <- prcomp(TrunkMid[c(1:7)], center = TRUE, scale = FALSE) # PCA

#"posterior" (12) comparative verts

TrunkPost <- rbind(T8_comb, T12_comb)
TrunkPost.pca <- prcomp(TrunkPost[c(1:7)], center = TRUE, scale = FALSE) # PCA





# MAKE FOSSIL VERTEBRAE GROUPINGS BASED ON ESTIMATED POSITION IN VERTEBRAL COLUMN #---------------------------------

Fossilverts <- dplyr::filter(PoAtVerts, grepl('41229*', species)) # fossils only
row.names(Fossilverts) <- Fossilverts$species
Fossilverts <- subset(Fossilverts, select=-c(specimen_num))

#Anterior trunk fossils

TrunkFossilT1<- dplyr::select(Fossilverts, contains("a", ignore.case = FALSE), contains("species"))
TrunkFossilT1 <- TrunkFossilT1 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)

TrunkFossilT4<- dplyr::select(Fossilverts, contains("b", ignore.case = FALSE), contains("species"))
TrunkFossilT4 <- TrunkFossilT4 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)

TrunkFossilAnt <- rbind(TrunkFossilT1, TrunkFossilT4)
TrunkFossilAnt <- na.omit(TrunkFossilAnt) # remove rows with N/A's

#Mid trunk fossils

TrunkFossilT8<- dplyr::select(Fossilverts, contains("c", ignore.case = FALSE), contains("species"))
TrunkFossilT8 <- TrunkFossilT8 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
TrunkFossilT8 <- na.omit(TrunkFossilT8) # remove rows with N/A's

#Posterior trunk fossils

TrunkFossilT12<- dplyr::select(Fossilverts, contains("d", ignore.case = FALSE), contains("species"))
TrunkFossilT12 <- TrunkFossilT12 %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
TrunkFossilT12 <- na.omit(TrunkFossilT12) # remove rows with N/A's

#Sacral fossils

TrunkFossilSc <- dplyr::select(Fossilverts, num_range("X", 14:20), contains("species"))
TrunkFossilSc <- TrunkFossilSc %>% dplyr::rename(X14 = 1, X15=2, X16=3,X17=4,X18=5,X19=6, X20=7)
TrunkFossilSc <- na.omit(TrunkFossilSc) # remove rows with N/A's






## PRINCIPAL COMPONENT ANALYSES WITH FOSSILS ##---------------------------------

# Anterior trunk fossils

TrunkFossilAnt_PCA <- predict(TrunkAnt.pca, TrunkFossilAnt[,1:7])
TrunkFossilAnt_PC_scores <- as.data.frame(TrunkFossilAnt_PCA)

FossilAnt_PC_scores <-data.frame(TrunkFossilAnt_PCA, species = TrunkFossilAnt$species)

Antscores <-data.frame(TrunkAnt.pca$x, species = TrunkAnt$species)

All_AntPC_scores <- (rbind(Antscores, FossilAnt_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of the fossils

# PLOT #

fossilcolors <- grDevices::gray.colors(54, start = 0, end = 0)
speciesshapes <- c(rep(16,15), rep(18,50))

library(ggplot2)
library(ggforce)

percentage_ant <- round(TrunkAnt.pca$sdev^2 / sum(TrunkAnt.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_ant <- paste( colnames(Antscores), "(", paste( as.character(percentage_ant), "%", ")", sep="") )

Ant_plot<-ggplot(All_AntPC_scores,aes(x=PC1,y=PC2,color=species, shape = species)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_ant[1]) + ylab(percentage_ant[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkAnt$species), values=c(speciescolors, fossilcolors)) + 
  scale_shape_manual(values = c(speciesshapes), guide = 'none') + theme_classic() + ggtitle("Anterior Vertebrae") + theme(legend.position = "none")
Ant_plot




# Middle trunk fossils

TrunkFossilMID_PCA <- predict(TrunkMid.pca, TrunkFossilT8[,1:7])
TrunkFossilMID_PC_scores <- as.data.frame(TrunkFossilMID_PCA)

MIDscores <-data.frame(TrunkMid.pca$x, species= (TrunkMid$species))

FossilMID_PC_scores <-data.frame(TrunkFossilMID_PCA, species= TrunkFossilT8$species)

All_MIDPC_scores <- (rbind(MIDscores, FossilMID_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil

percentage_mid <- round(TrunkMid.pca$sdev^2 / sum(TrunkMid.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_mid <- paste( colnames(MIDscores), "(", paste( as.character(percentage_mid), "%", ")", sep="") )

# PLOT #
Mid_plot<-ggplot(All_MIDPC_scores,aes(x=PC1,y=PC2,color=species, shape = species)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_mid[1]) + ylab(percentage_mid[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkMid$species), values=c(speciescolors, fossilcolors)) + 
  scale_shape_manual(values = c(speciesshapes), guide = 'none') +theme_classic() + ggtitle("Middle Vertebrae")+ theme(legend.position = "none")
Mid_plot




# Posterior trunk fossils

TrunkFossilPOST_PCA <- predict(TrunkPost.pca, TrunkFossilT12[,1:7])
TrunkFossilPOST_PC_scores <- as.data.frame(TrunkFossilPOST_PCA)

POSTscores <-data.frame(TrunkPost.pca$x, species= (TrunkPost$species))

FossilPOST_PC_scores <-data.frame(TrunkFossilPOST_PCA, species= TrunkFossilT12$species)

All_POSTPC_scores <- (rbind(POSTscores, FossilPOST_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil

percentage_post <- round(TrunkPost.pca$sdev^2 / sum(TrunkPost.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_post <- paste( colnames(POSTscores), "(", paste( as.character(percentage_post), "%", ")", sep="") )

# PLOT #
Post_plot<-ggplot(All_POSTPC_scores,aes(x=PC1,y=PC2,color=species, shape = species)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_post[1]) + ylab(percentage_post[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkPost$species), values=c(speciescolors, fossilcolors)) + 
  scale_shape_manual(values = c(speciesshapes), guide = 'none') + theme_classic() + ggtitle("Posterior Vertebrae")+ theme(legend.position = "none")
Post_plot




# Sacral fossils

FossilSC_PCA <- predict(Sc.pca, TrunkFossilSc[,1:7])
FossilSC_PC_scores <- as.data.frame(FossilSC_PCA)

SCscores <- data.frame(Sc.pca$x, species= Sc$species)


FossilSC_PC_scores <- data.frame(FossilSC_PCA, species= TrunkFossilSc$species)

All_SCPC_scores <- (rbind(SCscores, FossilSC_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil

percentage_Sc <- round(Sc.pca$sdev^2 / sum(Sc.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage_Sc <- paste( colnames(SCscores), "(", paste( as.character(percentage_Sc), "%", ")", sep="") )

# PLOT #
Sc_plot<-ggplot(All_SCPC_scores,aes(x=PC1,y=PC2,color=species, shape = species)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_Sc[1]) + ylab(percentage_Sc[2]) +
  scale_color_manual(name = "Species", breaks=levels(Sc$species), values=c(speciescolors, fossilcolors)) + 
  scale_shape_manual(values = c(speciesshapes), guide = 'none') +theme_classic() + ggtitle("Sacral Vertebrae")+ theme(legend.position = "none")
Sc_plot



# PLOT ALL TOGETHER #
par(oma=c(0,0,2,0))
gridExtra::grid.arrange(Ant_plot, Mid_plot ,Post_plot,Sc_plot,nrow = 2)






# Assess sample size per species---------------------------------
library(tidyverse)

PoAtVerts_wofossil %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

# *removed A. subsalsum and A. ordinarium* ---------------------------------

#Anterior vertebrae

TrunkAnt_sub <- dplyr::filter(TrunkAnt, !grepl('A.subsalsum|A.ordinarium', species))
TrunkAnt_sub$species <- droplevels(TrunkAnt_sub$species)

#Mid vertebrae
TrunkMid_sub <- dplyr::filter(TrunkMid, !grepl('A.subsalsum|A.ordinarium', species))
TrunkMid_sub$species <- droplevels(TrunkMid_sub$species)


#Posterior vertebrae
TrunkPost_sub <- dplyr::filter(TrunkPost, !grepl('A.subsalsum|A.ordinarium', species))
TrunkPost_sub$species <- droplevels(TrunkPost_sub$species)


# Sacral vertebrae
TrunkSc_sub <- dplyr::filter(Sc, !grepl('A.subsalsum|A.ordinarium', species))
TrunkSc_sub$species <- droplevels(TrunkSc_sub$species)





### CONSERVATIVE K NEAREST NEIGHBOR   ###:Non-parametric---------------------------------
library(caret)

#make KNN model using LOOCV to find optimal k

# SACRAL VERTEBRAE #


# LOOCV WITH REPLICATION
library(foreach)
library(doParallel)
ncore <- detectCores()
registerDoParallel(cores=ncore)

set.seed(123)
runs <- 1
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(species ~ X14 + X15 + X16 + X17 + X18 + X19 + X20,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkSc_sub)$results
#   }
# }) #repeated KNN model using LOOCV to find optimal k

fishSc <- map_dfr(fishSc,`[`, c("k", "Accuracy", "Kappa"))
kfishSc <- fishSc %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishSc <- kfishSc[1,] # k with highest accuracy


set.seed(123)
predicted.classes_sc <- train(species ~ X14 + X15 + X16 + X17 + X18 + X19 + X20,
            method     = "knn",
            tuneGrid   = expand.grid(k = 3),
            trControl  = trainControl(method  = "LOOCV"),
            metric     = "Accuracy",
            data       = TrunkSc_sub)$pred # predict class based on KNN model
mean(predicted.classes_sc$pred == predicted.classes_sc$obs) #overall accuracy

accKNNsc <- table(predicted.classes_sc$obs,predicted.classes_sc$pred)
accKNNsc
# t <- diag(prop.table(accKNNsc, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "Sacral species KNNAC", sep = ",", quote = FALSE, row.names = T)




# POSTERIOR VERTEBRAE #
set.seed(123)
runs <- 100
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkPost_sub)$results
#   }
# })

fishpost <- map_dfr(fishpost,`[`, c("k", "Accuracy", "Kappa"))
kfishpost <- fishpost %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishpost <- kfishpost[1,]

set.seed(123)
predicted.classes_post <- train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                                          method     = "knn",
                                          tuneGrid   = expand.grid(k = 1),
                                          trControl  = trainControl(method  = "LOOCV"),
                                          metric     = "Accuracy",
                                          data       = TrunkPost_sub)$pred # predict class based on KNN model
mean(predicted.classes_post$pred == predicted.classes_post$obs) #overall accuracy

accKNNpost <- table(predicted.classes_post$obs,predicted.classes_post$pred)
accKNNpost
# t <- diag(prop.table(accKNNpost, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "Posterior trunk linear species KNNAC", sep = ",", quote = FALSE, row.names = T)








# MIDDLE VERTEBRAE #
set.seed(123)
runs <- 100
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkMid_sub)$results
#   }
# })

fishmid <- map_dfr(fishmid,`[`, c("k", "Accuracy", "Kappa"))
kfishmid <- fishmid %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishmid <- kfishmid[1,]

set.seed(123)
predicted.classes_mid <- train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                                         method     = "knn",
                                         tuneGrid   = expand.grid(k = 1),
                                         trControl  = trainControl(method  = "LOOCV"),
                                         metric     = "Accuracy",
                                         data       = TrunkMid_sub)$pred # predict class based on KNN model
mean(predicted.classes_mid$pred == predicted.classes_mid$obs) #overall accuracy

accKNNmid <- table(predicted.classes_mid$obs,predicted.classes_mid$pred)
accKNNmid
# t <- diag(prop.table(accKNNmid, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "Middle trunk linear species KNNAC", sep = ",", quote = FALSE, row.names = T)




# ANTERIOR VERTEBRAE #
set.seed(123)
runs <- 100
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkAnt_sub)$results
#   }
# })

fishant <- map_dfr(fishant,`[`, c("k", "Accuracy", "Kappa"))
kfishant <- fishant %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishant <- kfishant[1,]

set.seed(123)
predicted.classes_ant <- train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                                         method     = "knn",
                                         tuneGrid   = expand.grid(k = 3),
                                         trControl  = trainControl(method  = "LOOCV"),
                                         metric     = "Accuracy",
                                         data       = TrunkAnt_sub)$pred # predict class based on KNN model
mean(predicted.classes_ant$pred == predicted.classes_ant$obs) #overall accuracy

accKNNant <- table(predicted.classes_ant$obs,predicted.classes_ant$pred)
accKNNant
# t <- diag(prop.table(accKNNant, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "Anterior verts linear species KNNAC", sep = ",", quote = FALSE, row.names = T)



## KNN FOSSIL CLASSIFICATIONS ##---------------------------------
library(class)

# ANTERIOR FOSSILS #

KnnAntPrediction <- knn(TrunkAnt_sub[,1:7], TrunkFossilAnt[,1:7],
                           TrunkAnt_sub$species, k= 3, prob=TRUE)
KnnAntPrediction
# t <- data.frame(fossil = TrunkFossilAnt$species, class = KnnAntPrediction)
# 
# write.table(t, file = "KnnAntPrediction fossils.txt", sep = ",", quote = FALSE, row.names = T)

# t <- cbind(as.character(TrunkFossilAnt$species), as.character(KnnAntPrediction_k7), as.character(attr(KnnAntPrediction_k7, 'prob')))
# t <- as.data.frame(t)
# t$V3<- as.numeric(t$V3)
# t$V3 <- round(t$V3, 2)
# write.table(t, file = "Anterior fossils.txt", sep = ",", quote = FALSE, row.names = T)


# MIDDLE FOSSILS #
KnnMidPrediction <- knn(TrunkMid_sub[,1:7], TrunkFossilT8[,1:7],
                           TrunkMid_sub$species, k=1, prob=TRUE)
KnnMidPrediction
# t <- data.frame(fossil = TrunkFossilT8$species, class = KnnMidPrediction)
# 
# write.table(t, file = "KnnMidPrediction fossils.txt", sep = ",", quote = FALSE, row.names = T)


# POSTERIOR FOSSILS #
KnnPostPrediction <- knn(TrunkPost_sub[,1:7], TrunkFossilT12[,1:7],
                            TrunkPost_sub$species, k=1, prob=TRUE)
KnnPostPrediction
# t <- data.frame(fossil = TrunkFossilT12$species, class = KnnPostPrediction)
# 
# write.table(t, file = "KnnPostPrediction fossils.txt", sep = ",", quote = FALSE, row.names = T)


# SACRAL VERTEBRAE FOSSILS #
KnnScPrediction <- knn(TrunkSc_sub[,1:7], TrunkFossilSc[,1:7],
                            TrunkSc_sub$species, k=3, prob=TRUE)
KnnScPrediction
# t <- data.frame(fossil = TrunkFossilSc$species, class = KnnScPrediction)
# 
# write.table(t, file = "KnnScPrediction fossils.txt", sep = ",", quote = FALSE, row.names = T)






### CONSERVATIVE RANDOM FOREST CLASSIFICATION ###:Non-parametric---------------------------------
library(randomForest)

# ANTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_ant <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkAnt_sub, importance=FALSE)
print(Atlas.rf_ant)
rf_acc_ant <- Atlas.rf_ant$confusion
rf_acc_ant <- 1-rf_acc_ant[,14] # percent correct classification
rf_acc_ant

# t <- rf_acc_ant
# t <-round(t, digits = 2)
# write.table(t, file = "Anterior verts RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_ant$predicted == TrunkAnt_sub$species) #overall accuracy



# FOSSIL CLASSIFICATION #
y_pred_ant = predict(Atlas.rf_ant, newdata = TrunkFossilAnt[,1:7])
y_pred_ant
# write.table(y_pred_ant, file = "Anterior fossils RF", sep = ",", quote = FALSE, row.names = T)



# MIDDLE VERTEBRAE #
set.seed(123)
Atlas.rf_mid <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkMid_sub, importance=FALSE)
print(Atlas.rf_mid)
rf_acc_mid <- Atlas.rf_mid$confusion
rf_acc_mid <- 1-rf_acc_mid[,14] # percent correct classification
rf_acc_mid
# t <- rf_acc_mid
# t <-round(t, digits = 2)
# write.table(t, file = "Middle verts RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_mid$predicted == TrunkMid_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_mid = predict(Atlas.rf_mid, newdata = TrunkFossilT8[,1:7])
y_pred_mid
# write.table(y_pred_mid, file = "Middle fossils RF", sep = ",", quote = FALSE, row.names = T)



# POSTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_post <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkPost_sub, importance=FALSE)
print(Atlas.rf_post)
rf_acc_post <- Atlas.rf_post$confusion
rf_acc_post <- 1-rf_acc_post[,14] # percent correct classification
rf_acc_post
# t <- rf_acc_post
# t <-round(t, digits = 2)
# write.table(t, file = "Posterior verts RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_post$predicted == TrunkPost_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_post = predict(Atlas.rf_post, newdata = TrunkFossilT12[,1:7])
y_pred_post
# write.table(y_pred_post, file = "Posterior fossils RF", sep = ",", quote = FALSE, row.names = T)



# SACRAL VERTEBRAE #
set.seed(123)
Atlas.rf_sc <- randomForest(species ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data=TrunkSc_sub, importance=FALSE)
print(Atlas.rf_sc)
rf_acc_sc <- Atlas.rf_sc$confusion
rf_acc_sc <- 1-rf_acc_sc[,14] # percent correct classification
rf_acc_sc
t <- rf_acc_sc
t <-round(t, digits = 2)
# write.table(t, file = "Sacral verts RFAC species", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_sc$predicted == TrunkSc_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_sc = predict(Atlas.rf_sc, newdata = TrunkFossilSc[,1:7])
y_pred_sc
# write.table(y_pred_sc, file = "Sacral fossils RF", sep = ",", quote = FALSE, row.names = T)





# ### MAHALAHOBIS DISTANCES AND NEIGHBOR JOINING ###---------------------------------
# library(HDMD)
# library("ape")
# library(phytools)
# 
# Mah_Dis <- function(pairwiseMah) {
#   names = rownames(pairwiseMah$means) #capture labels
#   mahala = sqrt(pairwiseMah$distance) #mahalanobis distance
#   rownames(mahala) = names #set rownames in the dissimilarity matrix
#   colnames(mahala) = names #set colnames in the dissimilarity matrix
#   return(mahala <- as.dist(mahala)) #this is the mahalanobis dissimilarity matrix 
# } # return mahalanobis dissimilarity matrix
# 
# 
# # ANTERIOR VERTEBRAE #
# ANTTotal <- rbind(TrunkAnt[,1:8], TrunkFossilAnt)
# Mahala1ANT = pairwise.mahalanobis(ANTTotal[,1:7], ANTTotal$species, digits = 3)
# Mah_DisAnt <- Mah_Dis(Mahala1ANT)
# 
# trANT <- nj(Mah_DisAnt) #neighbor joining
# 
# plot((as.phylo(trANT)),type="unrooted",cex=0.6,
#      use.edge.length=TRUE,lab4ut="axial",
#      no.margin=TRUE)
# 
# # MIDDLE VERTEBRAE #
# MIDTotal <- rbind(TrunkMid[,1:8], TrunkFossilT8)
# Mahala1MID = pairwise.mahalanobis(MIDTotal[,1:7], MIDTotal$species, digits = 3)
# Mah_DisMID <- Mah_Dis(Mahala1MID)
# 
# trMID <- nj(Mah_DisMID) #neighbor joining
# 
# plot((as.phylo(trMID)),type="unrooted",cex=0.6,
#      use.edge.length=TRUE,lab4ut="axial",
#      no.margin=TRUE)
# 
# 
# # POSTERIOR VERTEBRAE #
# POSTTotal <- rbind(TrunkPost[,1:8], TrunkFossilT12)
# Mahala1POST = pairwise.mahalanobis(POSTTotal[,1:7], POSTTotal$species, digits = 3)
# Mah_DisPOST <- Mah_Dis(Mahala1POST)
# 
# trPOST <- nj(Mah_DisPOST) #neighbor joining
# 
# plot((as.phylo(trPOST)),type="unrooted",cex=0.6,
#      use.edge.length=TRUE,lab4ut="axial",
#      no.margin=TRUE)# MIDDLE VERTEBRAE #
# 
# 
# # SACRAL VERTEBRAE #
# SCTotal <- rbind(Sc[,1:8], TrunkFossilSc)
# Mahala1SC = pairwise.mahalanobis(SCTotal[,1:7], SCTotal$species, digits = 3)
# Mah_DisSC <- Mah_Dis(Mahala1SC)
# 
# trSC <- nj(Mah_DisSC) #neighbor joining
# 
# plot((as.phylo(trSC)),type="unrooted",cex=0.6,
#      use.edge.length=TRUE,lab4ut="axial",
#      no.margin=TRUE)# MIDDLE VERTEBRAE #





# AMBYSTOMA CLADE CLASSIFICATION #---------------------------------
TrunkAnt_sub$clades <- dplyr::recode(TrunkAnt_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
TrunkMid_sub$clades <- dplyr::recode(TrunkMid_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
TrunkPost_sub$clades <- dplyr::recode(TrunkPost_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
TrunkSc_sub$clades <- dplyr::recode(TrunkSc_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")


### K NEAREST NEIGHBOR CLADES ###:Non-parametric---------------------------------
library(caret)

#make KNN model using LOOCV to find optimal k

# SACRAL VERTEBRAE #
library(caret)
# LOOCV WITH REPLICATION
library(foreach)
library(doParallel)
ncore <- detectCores()
registerDoParallel(cores=ncore)

set.seed(123)
runs <- 1
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(clades ~ X14 + X15 + X16 + X17 + X18 + X19 + X20,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkSc_sub)$results
#   }
# })

fishScclade <- map_dfr(fishScclade,`[`, c("k", "Accuracy", "Kappa"))
kfishScclade <- fishScclade %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishScclade <- kfishScclade[1,]


set.seed(123)
predicted.classes_sc_clade <- train(clades ~ X14 + X15 + X16 + X17 + X18 + X19 + X20,
                                              method     = "knn",
                                              tuneGrid   = expand.grid(k = 3),
                                              trControl  = trainControl(method  = "LOOCV"),
                                              metric     = "Accuracy",
                                              data       = TrunkSc_sub)$pred # predict class based on KNN model
mean(predicted.classes_sc_clade$pred == predicted.classes_sc_clade$obs) #overall accuracy

accKNN_sc_clade <- table(predicted.classes_sc_clade$obs,predicted.classes_sc_clade$pred)
accKNN_sc_clade
# t <- diag(prop.table(accKNN_sc_clade, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "Sacral clade KNNAC", sep = ",", quote = FALSE, row.names = T)




# POSTERIOR VERTEBRAE #
set.seed(123)
runs <- 1
# system.time({
#   fishpostCLADE <- foreach(icount(runs)) %dopar% {
#     train(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkPost_sub)$results
#   }
# })

fishpostCLADE <- map_dfr(fishpostCLADE,`[`, c("k", "Accuracy", "Kappa"))
kfishpostCLADE <- fishpostCLADE %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishpostCLADE <- kfishpostCLADE[1,]

set.seed(123)
predicted.classes_post_clade <- train(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                                                method     = "knn",
                                                tuneGrid   = expand.grid(k = 1),
                                                trControl  = trainControl(method  = "LOOCV"),
                                                metric     = "Accuracy",
                                                data       = TrunkPost_sub)$pred # predict class based on KNN model
mean(predicted.classes_post_clade$pred == predicted.classes_post_clade$obs) #overall accuracy

accKNNpostCLADE <- table(predicted.classes_post_clade$obs,predicted.classes_post_clade$pred)
accKNNpostCLADE
# t <- diag(prop.table(accKNNpostCLADE, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "posterior trunk clade KNNAC", sep = ",", quote = FALSE, row.names = T)






# MIDDLE VERTEBRAE #
set.seed(123)
runs <- 1
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkMid_sub)$results
#   }
# })

fishmid <- map_dfr(fishmid,`[`, c("k", "Accuracy", "Kappa"))
kfishmid <- fishmid %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishmid <- kfishmid[1,]

set.seed(123)
predicted.classes_mid <- train(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                                         method     = "knn",
                                         tuneGrid   = expand.grid(k = 1),
                                         trControl  = trainControl(method  = "LOOCV"),
                                         metric     = "Accuracy",
                                         data       = TrunkMid_sub)$pred # predict class based on KNN model
mean(predicted.classes_mid$pred == predicted.classes_mid$obs) #overall accuracy

accKNNmidCLADE <- table(predicted.classes_mid$obs,predicted.classes_mid$pred)
accKNNmidCLADE
# t <- diag(prop.table(accKNNmidCLADE, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "middelt trunk clade KNNAC", sep = ",", quote = FALSE, row.names = T)







# ANTERIOR VERTEBRAE #
set.seed(123)
runs <- 1
# system.time({
#   fishSc <- foreach(icount(runs)) %dopar% {
#     train(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
#           method     = "knn",
#           tuneGrid   = expand.grid(k = 1:17),
#           trControl  = trainControl(method  = "LOOCV"),
#           metric     = "Accuracy",
#           data       = TrunkAnt_sub)$results
#   }
# })

fishantCLADE <- map_dfr(fishantCLADE,`[`, c("k", "Accuracy", "Kappa"))
kfishantCLADE <- fishantCLADE %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kfishantCLADE <- kfishantCLADE[1,]

set.seed(123)
predicted.classes_ant_clade <- train(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a,
                                               method     = "knn",
                                               tuneGrid   = expand.grid(k = 6),
                                               trControl  = trainControl(method  = "LOOCV"),
                                               metric     = "Accuracy",
                                               data       = TrunkAnt_sub)$pred # predict class based on KNN model
mean(predicted.classes_ant_clade$pred == predicted.classes_ant_clade$obs) #overall accuracy

accKNN_ant_clade <- table(predicted.classes_ant_clade$obs,predicted.classes_ant_clade$pred)
accKNN_ant_clade
# t <- diag(prop.table(accKNN_ant_clade, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "anterior trunk clade KNNAC", sep = ",", quote = FALSE, row.names = T)








## KNN FOSSIL CLADE CLASSIFICATIONS ##---------------------------------
library(class)

# ANTERIOR FOSSILS #

KnnAntPrediction <- knn(TrunkAnt_sub[,1:7], TrunkFossilAnt[,1:7],
                        TrunkAnt_sub$clades, k=6 , prob=TRUE)
KnnAntPrediction
# t <- data.frame(fossil = TrunkFossilAnt$species, class = KnnAntPrediction)
# 
# write.table(t, file = "KnnAntPrediction clade fossils.txt", sep = ",", quote = FALSE, row.names = T)

# t <- cbind(as.character(TrunkFossilAnt$clades), as.character(KnnAntPrediction_k7), as.character(attr(KnnAntPrediction_k7, 'prob')))
# t <- as.data.frame(t)
# t$V3<- as.numeric(t$V3)
# t$V3 <- round(t$V3, 2)
# write.table(t, file = "Anterior fossils.txt", sep = ",", quote = FALSE, row.names = T)


# MIDDLE FOSSILS #
KnnMidPrediction <- knn(TrunkMid_sub[,1:7], TrunkFossilT8[,1:7],
                        TrunkMid_sub$clades, k=1, prob=TRUE)
KnnMidPrediction
# t <- data.frame(fossil = TrunkFossilT8$species, class = KnnMidPrediction)
# 
# write.table(t, file = "KnnMidPrediction clade fossils.txt", sep = ",", quote = FALSE, row.names = T)


# POSTERIOR FOSSILS #
KnnPostPrediction <- knn(TrunkPost_sub[,1:7], TrunkFossilT12[,1:7],
                         TrunkPost_sub$clades, k=1, prob=TRUE)
KnnPostPrediction
# t <- data.frame(fossil = TrunkFossilT12$species, class = KnnPostPrediction)
# 
# write.table(t, file = "KnnPostPrediction clade fossils.txt", sep = ",", quote = FALSE, row.names = T)


# SACRAL VERTEBRAE FOSSILS #
KnnSCPrediction <- knn(TrunkSc_sub[,1:7], TrunkFossilSc[,1:7],
                         TrunkSc_sub$clades, k=3, prob=TRUE)
KnnSCPrediction
# t <- data.frame(fossil = TrunkFossilSc$species, class = KnnSCPrediction)
# 
# write.table(t, file = "KnnSCPrediction clade fossils.txt", sep = ",", quote = FALSE, row.names = T)






### RANDOM FOREST CLADE CLASSIFICATION ###:Non-parametric---------------------------------
library(randomForest)

# ANTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_ant_clade <- randomForest(clades ~ X1a + X2a + X3a + X4a + X5a + X6a 
                                   # + X7a
                                   , data=TrunkAnt_sub, importance=FALSE)
print(Atlas.rf_ant_clade)
rf_acc_ant_clade <- Atlas.rf_ant_clade$confusion
rf_acc_ant_clade <- 1-rf_acc_ant_clade[,9] # percent correct classification
rf_acc_ant_clade
# t <- rf_acc_ant_clade
# t <-round(t, digits = 2)
# write.table(t, file = "Anterior verts RFAC clade", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_ant_clade$predicted == TrunkAnt_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_ant_clade = predict(Atlas.rf_ant_clade, newdata = TrunkFossilAnt[,1:7])
y_pred_ant_clade
# write.table(y_pred_ant_clade, file = "ant clade fossils RF.txt", sep = ",", quote = FALSE, row.names = T)



# MIDDLE VERTEBRAE #
set.seed(123)
Atlas.rf_mid_clade <- randomForest(clades ~ X1a + X2a + X3a + X4a + X5a + X6a 
                                   # + X7a
                                   , data=TrunkMid_sub, importance=FALSE)
print(Atlas.rf_mid_clade)
rf_acc_mid_clade <- Atlas.rf_mid_clade$confusion
rf_acc_mid_clade <- 1-rf_acc_mid_clade[,9] # percent correct classification
rf_acc_mid_clade
# t <- rf_acc_mid_clade
# t <-round(t, digits = 2)
# write.table(t, file = "Middle verts RFAC clade", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_mid_clade$predicted == TrunkMid_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_mid_clade = predict(Atlas.rf_mid_clade, newdata = TrunkFossilT8[,1:7])
y_pred_mid_clade
# write.table(y_pred_mid_clade, file = "middle clade fossils RF.txt", sep = ",", quote = FALSE, row.names = T)



# POSTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_post_clade <- randomForest(clades ~ X1a + X2a + X3a + X4a + X5a + X6a 
                                    # + X7a
                                    , data=TrunkPost_sub, importance=FALSE)
print(Atlas.rf_post_clade)
rf_acc_post_clade <- Atlas.rf_post_clade$confusion
rf_acc_post_clade <- 1-rf_acc_post_clade[,9] # percent correct classification
rf_acc_post_clade
# t <- rf_acc_post_clade
# t <-round(t, digits = 2)
# write.table(t, file = "Posterior verts RFAC clade", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_post_clade$predicted == TrunkPost_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_post_clade = predict(Atlas.rf_post_clade, newdata = TrunkFossilT12[,1:7])
y_pred_post_clade
# write.table(y_pred_post_clade, file = "post clade fossils RF.txt", sep = ",", quote = FALSE, row.names = T)



# SACRAL VERTEBRAE #
set.seed(123)
Atlas.rf_sc_clade <- randomForest(clades ~ X14 + X15 + X16 + X17 + X18 + X19 
                                  # + X20
                                  , data=TrunkSc_sub, importance=FALSE)
print(Atlas.rf_sc_clade)
rf_acc_sc_clade <- Atlas.rf_sc_clade$confusion
rf_acc_sc_clade <- 1-rf_acc_sc_clade[,9] # percent correct classification
rf_acc_sc_clade
# t <- rf_acc_sc_clade
# t <-round(t, digits = 2)
# write.table(t, file = "Sacral verts RFAC clade", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_sc_clade$predicted == TrunkSc_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_sc_clade = predict(Atlas.rf_sc_clade, newdata = TrunkFossilSc[,1:7])
y_pred_sc_clade
# write.table(y_pred_sc_clade, file = "sc clade fossils RF.txt", sep = ",", quote = FALSE, row.names = T)











### SEXUAL DIMORPHISM? ###---------------------------------

PoAtVerts_sex <-  Amb_linear_data[c(6:7, 14:48, 53:55)] #select only relevant Post atlantal measurements

# REMOVE FOSSILS #
PoAtVerts_sex <- dplyr::filter(PoAtVerts_sex, !grepl('41229*', species)) # remove fossils
row.names(PoAtVerts_sex) <- PoAtVerts_sex$specimen_num
# PoAtVerts_wofossil <- subset(PoAtVerts_wofossil, select=-c(specimen_num))
PoAtVerts_sex <- droplevels(PoAtVerts_sex)
PoAtVerts_sex$species <- factor(PoAtVerts_sex$species, levels = 
                                  c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                    "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels
library(tidyverse)
# SUBSET DATA FOR EACH VERTEBRA #
T1_sex <- dplyr::select(PoAtVerts_sex, contains("a", ignore.case = FALSE),  contains("species"), contains("sex")) %>% drop_na()

T4_sex <- dplyr::select(PoAtVerts_sex, contains("b", ignore.case = FALSE), contains("species"), contains("sex")) %>% drop_na()

T8_sex <- dplyr::select(PoAtVerts_sex, contains("c", ignore.case = FALSE), contains("species"), contains("sex")) %>% drop_na()
T8_sex <- subset(T8_sex, select=-c(specimen_num))

T8_extension_sex <- dplyr::select(PoAtVerts_sex, contains("Cen"), contains("species"), contains("sex")) %>% drop_na()

T12_sex<- dplyr::select(PoAtVerts_sex, contains("d", ignore.case = FALSE), contains("species"), contains("sex")) %>% drop_na()

Sc_sex <- dplyr::select(PoAtVerts_sex, num_range("X", 14:20), contains("species"), contains("sex")) %>% drop_na()
Sc.pca_sex <- prcomp(Sc_sex[c(1:7)], center = TRUE, scale = FALSE) # PCA





# MAKE REGIONAL VERTEBRAE GROUPING SUBDATASETS "ANTERIOR, MIDDLE, POSTERIOR" #

#"Anterior"(1,4) comparative verts

T1_comb_sex <- T1_sex %>%                          # Applying row_number function
  dplyr::mutate(Vert = "a")
T1_comb_sex <- T1_comb_sex %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)

T4_comb_sex <- T4_sex %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T4_comb_sex <- T4_comb_sex %>%                          # Applying row_number function
  dplyr::mutate(Vert = "b")

TrunkAnt_sex <- rbind(T1_comb_sex, T4_comb_sex)
TrunkAnt.pca_sex <- prcomp(TrunkAnt_sex[c(1:7)], center = TRUE, scale = FALSE) # PCA

#"Mid"(4,8,12) comparative verts
T8_comb_sex <- T8_sex %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T8_comb_sex <- T8_comb_sex %>%                          # Applying row_number function
  dplyr::mutate(Vert = "c")

T12_comb_sex <- T12_sex %>% dplyr::rename(X1a = 1, X2a=2, X3a=3,X4a=4,X5a=5,X6a=6, X7a=7)
T12_comb_sex <- T12_comb_sex %>%                          # Applying row_number function
  dplyr::mutate(Vert = "d")

TrunkMid_sex <- rbind(T4_comb_sex, T8_comb_sex, T12_comb_sex)
TrunkMid.pca_sex <- prcomp(TrunkMid_sex[c(1:7)], center = TRUE, scale = FALSE) # PCA

#"posterior" (12) comparative verts

TrunkPost_sex <- rbind(T8_comb_sex, T12_comb_sex)
TrunkPost.pca_sex <- prcomp(TrunkPost_sex[c(1:7)], center = TRUE, scale = FALSE) # PCA





library(ggplot2)
library(ggforce)

Antscores_sex <-data.frame(TrunkAnt.pca_sex$x, species = TrunkAnt_sex$species, sex = TrunkAnt_sex$Sex)

percentage_ant_sex <- round(TrunkAnt.pca_sex$sdev / sum(TrunkAnt.pca_sex$sdev) * 100, 2)# find percentage variance explained by PC's
percentage_ant_sex <- paste( colnames(Antscores_sex), "(", paste( as.character(percentage_ant_sex), "%", ")", sep="") )

speciescolors <- c("#666600", "#C9D42D" ,"#42CEBC" ,"#F2AB1F" ,"#864ED0" ,"#261ACE", "#086AFD", "#08FD6A", "#0C8B3F", "#E50CF5", "#FF5E00","#FF0000", "#FF6A6A", "#D5930F", "#9E1616")


Ant_plot_sex<-ggplot(Antscores_sex,aes(x=PC1,y=PC2,color=species, shape = sex)) + 
  #geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_ant_sex[1]) + ylab(percentage_ant_sex[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkAnt_sex$species), values=c(speciescolors)) + 
  theme_classic() + ggtitle("Anterior Vertebrae")
Ant_plot_sex




