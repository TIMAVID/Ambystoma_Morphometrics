## Load in data ##
library(readxl)
download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/Linear_data/Data/Amb_linear_data.xlsx?raw=true",   "Amb_linear_data.xlsx")
Amb_linear_data <- read_excel("Amb_linear_data.xlsx")

### Post atlantal vertebrae ###

# Tidy data #
library(dplyr)

PoAtVerts <-  Amb_linear_data[c(7:8, 15:49, 55:56)] #select only relevant Post atlantal measurements

PoAtVerts_wofossil <- dplyr::filter(PoAtVerts, !grepl('41229*', species)) # remove fossils
PoAtVerts_wofossil <- PoAtVerts_wofossil %>% dplyr::rename(M14 = 31, M15=32, M16=33,M17=34,M18=35,M19=36,M20=37) # remame columns for easier manipulation

PoAtVerts_wofossil_noNA <- na.omit(PoAtVerts_wofossil) # remove rows with N/A's

PoAtVerts_wofossil_noNA$species<-as.factor(PoAtVerts_wofossil_noNA$species)
PoAtVerts_wofossil_noNA$species <- factor(PoAtVerts_wofossil_noNA$species, levels = 
                                        c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                          "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")) # Reorder species
levels(PoAtVerts_wofossil_noNA$species)


T1 <- dplyr::select(PoAtVerts_wofossil_noNA, contains("a", ignore.case = FALSE), contains("specimen"), contains("species"))
T1 <- T1 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)

T4 <- dplyr::select(PoAtVerts_wofossil_noNA, contains("b", ignore.case = FALSE), contains("specimen"), contains("species"))
T4 <- T4 %>% dplyr::rename(M1b = 1, M2b=2, M3b=3,M4b=4,M5b=5,M6b=6, M7b=7)

T8 <- dplyr::select(PoAtVerts_wofossil_noNA, contains("c", ignore.case = FALSE), contains("specimen"), contains("species"))
T8 <- T8 %>% dplyr::rename(M1c = 1, M2c=2, M3c=3,M4c=4,M5c=5,M6c=6, M7c=7)

T8_extension <- dplyr::select(PoAtVerts_wofossil_noNA, contains("Cen"), contains("specimen"), contains("species"))

T12<- dplyr::select(PoAtVerts_wofossil_noNA, contains("d", ignore.case = FALSE), contains("specimen"), contains("species"))
T12 <- T12 %>% dplyr::rename(M1d = 1, M2d=2, M3d=3,M4d=4,M5d=5,M6d=6, M7d=7)

Sc <- dplyr::select(PoAtVerts_wofossil_noNA, contains("M", ignore.case = FALSE), contains("specimen"), contains("species"))


#"Anterior"(1,4) comparative verts

T1_comb <- T1 %>%                          # Applying row_number function
  dplyr::mutate(Vert = "a")

T4_comb <- T4 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)
T4_comb <- T4_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "b")

TrunkAnt <- rbind(T1_comb, T4_comb)
TrunkAnt.pca <- prcomp(TrunkAnt[c(1:7)], center = TRUE, scale = FALSE) # PCA

#"Mid"(4,8,12) comparative verts
T8_comb <- T8 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)
T8_comb <- T8_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "c")

T12_comb <- T12 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)
T12_comb <- T12_comb %>%                          # Applying row_number function
  dplyr::mutate(Vert = "d")

TrunkMid <- rbind(T4_comb, T8_comb, T12_comb)
TrunkMid.pca <- prcomp(TrunkMid[c(1:7)], center = TRUE, scale = FALSE) # PCA

#"posterior" (12) comparative verts

TrunkPost <- rbind(T8_comb, T12_comb)
TrunkPost.pca <- prcomp(TrunkPost[c(1:7)], center = TRUE, scale = FALSE) # PCA


#Anterior trunk fossils

TrunkFossilT1<- dplyr::select(PoAtVerts, contains("a", ignore.case = FALSE), contains("species"))
TrunkFossilT1 <- dplyr::filter(TrunkFossilT1, grepl('41229*', species)) # fossils
TrunkFossilT1 <- TrunkFossilT1 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)

TrunkFossilT4<- dplyr::select(PoAtVerts, contains("b", ignore.case = FALSE), contains("species"))
TrunkFossilT4 <- TrunkFossilT4 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)
TrunkFossilT4 <- dplyr::filter(TrunkFossilT4, grepl('41229*', species)) # fossils

TrunkFossilAnt <- rbind(TrunkFossilT1, TrunkFossilT4)
TrunkFossilAnt <- na.omit(TrunkFossilAnt) # remove rows with N/A's

#Mid trunk fossils

TrunkFossilT8<- dplyr::select(PoAtVerts, contains("c", ignore.case = FALSE), contains("species"))
TrunkFossilT8<- dplyr::select(TrunkFossilT8, -specimen_num)
TrunkFossilT8 <- TrunkFossilT8 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)
TrunkFossilT8 <- dplyr::filter(TrunkFossilT8, grepl('41229*', species)) # fossils
TrunkFossilT8 <- na.omit(TrunkFossilT8) # remove rows with N/A's

#Posterior trunk fossils

TrunkFossilT12<- dplyr::select(PoAtVerts, contains("d", ignore.case = FALSE), contains("species"))
TrunkFossilT12 <- TrunkFossilT12 %>% dplyr::rename(M1a = 1, M2a=2, M3a=3,M4a=4,M5a=5,M6a=6, M7a=7)
TrunkFossilT12 <- dplyr::filter(TrunkFossilT12, grepl('41229*', species)) # fossils
TrunkFossilT12 <- na.omit(TrunkFossilT12) # remove rows with N/A's

#Sacral fossils

TrunkFossilSc<- PoAtVerts[c(31:37,39)]
TrunkFossilSc <- TrunkFossilSc %>% dplyr::rename(M14 = 1, M15=2, M16=3,M17=4,M18=5,M19=6, M20=7)
TrunkFossilSc <- dplyr::filter(TrunkFossilSc, grepl('41229*', species)) # fossils
TrunkFossilSc <- na.omit(TrunkFossilSc) # remove rows with N/A's

## PCA's ##

T1.pca <- prcomp(T1[c(1:7)], center = TRUE, scale = FALSE) # PCA
T4.pca <- prcomp(T4[c(1:7)], center = TRUE, scale = FALSE) # PCA
T8.pca <- prcomp(T8[c(1:7)], center = TRUE, scale = TRUE) # PCA
T12.pca <- prcomp(T1[c(1:7)], center = TRUE, scale = FALSE) # PCA
Sc.pca <- prcomp(Sc[c(1:7)], center = TRUE, scale = FALSE) # PCA

library(ggbiplot)

par(mfrow=c(3,2))

ggbiplot(T1.pca, ellipse=FALSE, groups=T1$species)
ggbiplot(T4.pca, ellipse=FALSE, groups=T4$species)
ggbiplot(T8.pca, ellipse=FALSE, groups=T8$species)
ggbiplot(T12.pca, ellipse=FALSE, groups=T12$species)
ggbiplot(Sc.pca, ellipse=FALSE, groups=Sc$species)


library(ggfortify)
library(ggplot2)
library(gridExtra)
T1pcaPlot<- autoplot(T1.pca, data = T1, colour = 'species', frame = TRUE, label = TRUE) + ggtitle("T1")
T4pcaPlot<-  autoplot(T4.pca, data = T4, colour = 'species', frame = TRUE, label = TRUE) + ggtitle("T4")
T8pcaPlot<-  autoplot(T8.pca, data = T8, colour = 'species', frame = TRUE, label = TRUE) + ggtitle("T8")
T12pcaPlot<- autoplot(T12.pca, data = T12, colour = 'species', frame = TRUE, label = TRUE) + ggtitle("T12")
ScpcaPlot<-  autoplot(Sc.pca, data = Sc, colour = 'species', frame = TRUE, label = TRUE) + ggtitle("Sc")

grid.arrange(T1pcaPlot, T4pcaPlot,T8pcaPlot ,T12pcaPlot,ScpcaPlot,nrow = 3)



## PCA w/ fossils

#Anterior fossils
TrunkFossilAnt_PCA <- predict(TrunkAnt.pca, TrunkFossilAnt[,1:7])
TrunkFossilAnt_PC_scores <- as.data.frame(TrunkFossilAnt_PCA)

Antscores <-as.data.frame(TrunkAnt.pca$x)

Antscores <-cbind((Antscores), species= (TrunkAnt$species))

FossilAnt_PC_scores <-as.data.frame(TrunkFossilAnt_PCA)
FossilAnt_PC_scores <- cbind(FossilAnt_PC_scores, species= TrunkFossilAnt$species)

All_AntPC_scores <- (rbind(Antscores, FossilAnt_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil


AntFossils <- as.character(TrunkFossilAnt$species)

species <- c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
             "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")

library(ggalt)


ANTpcaplot <- ggplot(data = All_AntPC_scores, mapping = aes(x = PC1, y = PC2,
                                                      col = species, label = species)) # creates the initial plot with datapoints color-coded and unique symbols by each species
ANTpcaplot <- ANTpcaplot + geom_encircle(expand=0, size = 2, data = All_AntPC_scores[!All_AntPC_scores$species %in% AntFossils,])+ theme_bw()
ANTpcaplot <- ANTpcaplot + 
  geom_text(aes(PC1, PC2, label = species), nudge_y = .003,
            check_overlap = FALSE, data = All_AntPC_scores[All_AntPC_scores$species %in% AntFossils,])+ geom_point(data = All_AntPC_scores[All_AntPC_scores$species %in% AntFossils,])

ANTpcaplot <- ANTpcaplot + scale_color_manual(breaks = c(species),
                                        values=c("black", "black","black", "black","black", "black","black", "black","black", "black","black", "black","black", "black","black", "#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                                 "#CDDADD", "#DC5956", "#E363BB"))
ANTpcaplot


#Mid trunk fossils

TrunkFossilMID_PCA <- predict(TrunkMid.pca, TrunkFossilT8[,1:7])
TrunkFossilMID_PC_scores <- as.data.frame(TrunkFossilMID_PCA)

MIDscores <-as.data.frame(TrunkMid.pca$x)

MIDscores <-cbind((MIDscores), species= (TrunkMid$species))

FossilMID_PC_scores <-data.frame(TrunkFossilMID_PCA)
FossilMID_PC_scores <- cbind(FossilMID_PC_scores, species= TrunkFossilT8$species)

All_MIDPC_scores <- (rbind(MIDscores, FossilMID_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil


MIDFossils <- as.character(TrunkFossilT8$species)


MIDpcaplot <- ggplot(data = All_MIDPC_scores, mapping = aes(x = PC1, y = PC2,
                                                         col = species, label = species)) # creates the initial plot with datapoints color-coded and unique symbols by each species
MIDpcaplot <- MIDpcaplot + geom_encircle(expand=0, size = 2, data = All_MIDPC_scores[!All_MIDPC_scores$species %in% MIDFossils,])+ theme_bw()
MIDpcaplot <- MIDpcaplot + 
  geom_text(aes(PC1, PC2, label = species), nudge_y = .003,
            check_overlap = TRUE, data = All_MIDPC_scores[All_MIDPC_scores$species %in% MIDFossils,])+ geom_point(data = All_MIDPC_scores)

library(grDevices)
n <- grDevices::gray.colors(41, start = 0.1, end = 0.4)

MIDpcaplot <- MIDpcaplot + scale_color_manual(breaks = c(species),
                                        values=c(n, "#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                                 "#CDDADD", "#DC5956", "#E363BB"))
MIDpcaplot


#Posterior trunk fossils

TrunkFossilPOST_PCA <- predict(TrunkPost.pca, TrunkFossilT12[,1:7])
TrunkFossilPOST_PC_scores <- as.data.frame(TrunkFossilPOST_PCA)

POSTscores <-as.data.frame(TrunkPost.pca$x)

POSTscores <-cbind((POSTscores), species= (TrunkPost$species))

FossilPOST_PC_scores <-data.frame(TrunkFossilPOST_PCA)
FossilPOST_PC_scores <- cbind(FossilPOST_PC_scores, species= TrunkFossilT12$species)

All_POSTPC_scores <- (rbind(POSTscores, FossilPOST_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil


POSTFossils <- as.character(TrunkFossilT12$species)


POSTpcaplot <- ggplot(data = All_POSTPC_scores, mapping = aes(x = PC1, y = PC2,
                                                            col = species, label = species)) # creates the initial plot with datapoints color-coded and unique symbols by each species
POSTpcaplot <- POSTpcaplot + geom_encircle(expand=0, size = 2, data = All_POSTPC_scores[!All_POSTPC_scores$species %in% POSTFossils,])+ theme_bw()
POSTpcaplot <- POSTpcaplot + 
  geom_text(aes(PC1, PC2, label = species), nudge_y = .003,
            check_overlap = FALSE, data = All_POSTPC_scores[All_POSTPC_scores$species %in% POSTFossils,])+ geom_point(data = All_POSTPC_scores)

library(grDevices)
n <- grDevices::gray.colors(3, start = 0.1, end = 0.4)

POSTpcaplot <- POSTpcaplot + scale_color_manual(breaks = c(species),
                                              values=c(n, "#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                                       "#CDDADD", "#DC5956", "#E363BB"))
POSTpcaplot

# Sacral fossils

FossilSC_PCA <- predict(Sc.pca, TrunkFossilSc[,1:7])
FossilSC_PC_scores <- as.data.frame(FossilSC_PCA)

SCscores <-as.data.frame(Sc.pca$x)

SCscores <-cbind((SCscores), species= (Sc$species))

FossilSC_PC_scores <- data.frame(FossilSC_PCA)
FossilSC_PC_scores <- cbind(FossilSC_PC_scores, species= TrunkFossilSc$species)

All_SCPC_scores <- (rbind(SCscores, FossilSC_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil


SCFossils <- as.character(TrunkFossilSc$species)


SCpcaplot <- ggplot(data = All_SCPC_scores, mapping = aes(x = PC1, y = PC2,
                                                              col = species, label = species)) # creates the initial plot with datapoints color-coded and unique symbols by each species
SCpcaplot <- SCpcaplot + geom_encircle(expand=0, size = 2, data = All_SCPC_scores[!All_SCPC_scores$species %in% SCFossils,])+ theme_bw()
SCpcaplot <- SCpcaplot + 
  geom_text(aes(PC1, PC2, label = species), nudge_y = .003,
            check_overlap = FALSE, data = All_SCPC_scores[All_SCPC_scores$species %in% SCFossils,])+ geom_point(data = All_SCPC_scores)

library(grDevices)
n <- grDevices::gray.colors(7, start = 0.1, end = 0.4)

SCpcaplot <- SCpcaplot + scale_color_manual(breaks = c(species),
                                                values=c(n, "#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                                         "#CDDADD", "#DC5956", "#E363BB"))
SCpcaplot


# Posterior extension plot

library(EnvStats)
T8_extension

T8_extension<-mutate(T8_extension ,ratio = Cen_to_NeuAr/Cen_to_PoZy, .before = Cen_to_PoZy) #add new column for ratio

library(tidyr)
data_long <- gather(T8_extension, Type , Measurement, Cen_to_NeuAr:Cen_to_PoZy, factor_key=TRUE) #convert to long format

library(ggplot2)
s <- ggplot(data_long, aes(species, Measurement, fill = Type)) + geom_boxplot(position = "dodge")
s


# Assess sample size per species
library(tidyverse)

PoAtVerts_wofossil_noNA %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

T1_sub <- dplyr::filter(T1, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
T1_sub$species <- droplevels(T1_sub$species)
T1_sub <- tibble::column_to_rownames(T1_sub, var = "specimen_num")

T4_sub <- dplyr::filter(T4, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
T4_sub$species <- droplevels(T4_sub$species)
T4_sub <- tibble::column_to_rownames(T4_sub, var = "specimen_num")

T8_sub <- dplyr::filter(T8, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
T8_sub$species <- droplevels(T8_sub$species)
T8_sub <- tibble::column_to_rownames(T8_sub, var = "specimen_num")

T12_sub <- dplyr::filter(T12, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
T12_sub$species <- droplevels(T12_sub$species)
T12_sub <- tibble::column_to_rownames(T12_sub, var = "specimen_num")

Sc_sub <- dplyr::filter(Sc, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
Sc_sub$species <- droplevels(Sc_sub$species)
Sc_sub <- tibble::column_to_rownames(Sc_sub, var = "specimen_num")


### Random Forest ###:Non-parametric

library(randomForest)

T1.rf <- randomForest(species ~ M1a + M2a+ M3a+M4a+M5a+M6a+M7a, data=T1_sub, importance=TRUE,proximity=TRUE)
T4.rf <- randomForest(species ~ M1b + M2b+ M3b+M4b+M5b+M6b+M7b, data=T4_sub, importance=TRUE,proximity=TRUE)
T8.rf <- randomForest(species ~ M1c + M2c+ M3c+M4c+M5c+M6c+M7c, data=T8_sub, importance=TRUE,proximity=TRUE)
T12.rf <- randomForest(species ~ M1d + M2d+ M3d+M4d+M5d+M6d+M7d, data=T12_sub, importance=TRUE,proximity=TRUE)
Sc.rf <- randomForest(species ~ M14 + M15+ M16+M17+M18+M19+M20, data=Sc_sub, importance=TRUE,proximity=TRUE)

RFacc <- function(rfmodel){
  i <- rfmodel$confusion
  i <- 1-i[,11] # percent correct classification
  print(i)
  print(rfmodel)
}

RFacc(T1.rf)
RFacc(T4.rf)
RFacc(T8.rf)
RFacc(T12.rf)
RFacc(Sc.rf)

# Look at variable importance
round(importance(T1.rf), 2)
varImpPlot(T1.rf)


### Predict fossils ###

#Anterior fossils

TrunkAnt_sub <- dplyr::filter(TrunkAnt, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
TrunkAnt_sub$species <- droplevels(TrunkAnt_sub$species)
Ant.rf <- randomForest(species ~ M1a + M2a+ M3a+M4a+M5a+M6a+M7a, data=TrunkAnt_sub, importance=TRUE,proximity=TRUE)
AntFossil_pred = predict(Ant.rf, newdata = TrunkFossilAnt[,1:7])
AntFossil_pred

#Mid fossils
TrunkMid_sub <- dplyr::filter(TrunkMid, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
TrunkMid_sub$species <- droplevels(TrunkMid_sub$species)
Mid.rf <- randomForest(species ~ M1a + M2a+ M3a+M4a+M5a+M6a+M7a, data=TrunkMid, importance=TRUE,proximity=TRUE)
MidFossil_pred = predict(Mid.rf, newdata = TrunkFossilT8[,1:7])
MidFossil_pred

#Posterior fossils
TrunkPost_sub <- dplyr::filter(TrunkPost, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
TrunkPost_sub$species <- droplevels(TrunkPost_sub$species)
Post.rf <- randomForest(species ~ M1a + M2a+ M3a+M4a+M5a+M6a+M7a, data=TrunkPost, importance=TRUE,proximity=TRUE)
PostFossil_pred = predict(Post.rf, newdata = TrunkFossilT12[,1:7])
PostFossil_pred

# Sacral vert fossils
TrunkSc_sub <- dplyr::filter(Sc, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
TrunkSc_sub$species <- droplevels(TrunkSc_sub$species)
Sc.rf <- randomForest(species ~ M14 + M15+ M16+M17+M18+M19+M20, data=Sc, importance=TRUE,proximity=TRUE)
ScFossil_pred = predict(Sc.rf, newdata = TrunkFossilSc[,1:7])
ScFossil_pred


### K Nearest neighbor ###:Non-parametric

library(caret)

#make KNN model using LOOCV to find optimal k

set.seed(123)

KNNmodel <- train(
  species ~M14 + M15+ M16+M17+M18+M19+M20, data = TrunkSc_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(TrunkMid_sub[,1:7]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == TrunkMid_sub$species) #overall accuracy

# assess accuracy per species
accKNN <- table(TrunkMid_sub$species,predicted.classes)
accKNN
diag(prop.table(accKNN, 1))


## Fossil predictions ##

# Ant fossils #
library(class)
KnnAntPrediction_k9 <- knn(TrunkAnt_sub[,1:7], TrunkFossilAnt[,1:7],
                           TrunkAnt_sub$species, k=9, prob=TRUE)
KnnAntPrediction_k9

# Mid fossils #
KnnMidPrediction_k9 <- knn(TrunkMid_sub[,1:7], TrunkFossilT8[,1:7],
                           TrunkMid_sub$species, k=9, prob=TRUE)
KnnMidPrediction_k9

# Posterior fossils #
KnnPostPrediction_k9 <- knn(TrunkPost_sub[,1:7], TrunkFossilT12[,1:7],
                            TrunkPost_sub$species, k=9, prob=TRUE)
KnnPostPrediction_k9

# Sacral fossils #
KnnPostPrediction_k9 <- knn(TrunkSc_sub[,1:7], TrunkFossilSc[,1:7],
                            TrunkSc_sub$species, k=5, prob=TRUE)
KnnPostPrediction_k9


### Model selection (multinominal regression) 

library(glmulti)

library(nnet)

multinom.glmulti <- function(formula, data, ...)
  multinom(formula, data, ...)

res <- glmulti(species ~ .,level=1, data=T8_sub, report = FALSE, plotty = FALSE,fitfunction=multinom.glmulti, method = "h", crit="aic", confsetsize=16)

print(res)

plot(res)

top <- weightable(res)
top


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
Amb_species<-unique(T8_sub$species)
tips<-tree$tip.label
ii<-sapply(Amb_species,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
plotTree(tree,ftype="i")
#Tree did not include A.mavortium so I lumped that species with A.tigrinum
T8_sub_tig<- T8_sub
T8_sub_tig$species<-gsub("A.mavortium", "A.tigrinum", T8_sub_tig$species, fixed = TRUE)
T8_sub_tig$species<-as.factor(T8_sub_tig$species)
T8_sub_tig$species <- factor(T8_sub_tig$species, levels = 
                                             c("A.gracile", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum",
                                               "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium")) # Reorder species

#Preformed a group PCA
library(Morpho)
library(geomorph)
gpca <- groupPCA(T8_sub_tig[,1:7], T8_sub_tig$species, rounds=0)
plot(gpca$groupmeans)

#Performed a Phylogenetic PCA based on group means
phylo.PCA <- gm.prcomp(gpca$groupmeans, phy = tree, align.to.phy = FALSE)
summary(phylo.PCA)

A_species<-attributes(gpca$groupmeans) #access attributes names
A_species<-(A_species$dimnames[[1]])
A_species<-as.factor(A_species)
A_species <- factor(A_species, levels = 
                      c("A.gracile", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum",
                        "A.mabeei","A.texanum","A.annulatum","A.tigrinum")) # Reorder species

#Plot phylogenetic PCA
plot(phylo.PCA, phylo = TRUE, main = "phylo PCA", col=A_species)

#Test for phylogenetic signal, uses Blomberg’s K to test for strength and significance of phylogenetic signal.
physignal(gpca$groupmeans, tree, print.progress = F, iter = 999)

#Phylogenetic generalized least squares

avg_gdf<-geomorph.data.frame(coords=gpca$groupmeans, species=A_species) #make new geomorph dataframe with group mean coords

pgls<-procD.pgls(coords~species, phy=tree, data=avg_gdf, print.progress = F, iter = 999) #Phylogenetic generalized least squares
pgls$aov.table


#Compare evolutionary rates in different portions of the tree based on brownian motion

names(A_species) <- levels(A_species)

rate.comp<-compare.evol.rates(avg_gdf$coords, tree, gp=A_species, method = c("permutation"), iter = 999, print.progress = F)

plot(rate.comp)
rate.comp$sigma.d.gp
rate.comp$pairwise.pvalue

