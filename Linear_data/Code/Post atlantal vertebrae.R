## Load in data ##
library(readxl)
download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/Linear_data/Data/Amb_linear_data.xlsx?raw=true",   "Amb_linear_data.xlsx")
Amb_linear_data <- read_excel("Amb_linear_data.xlsx")

### Post atlantal vertebrae ###

# Tidy data #
library(dplyr)

PoAtVerts <-  Amb_linear_data[c(7:8, 15:49, 55:56)] #select only relevant Post atlantal measurements

PoAtVerts_wofossil <- dplyr::filter(PoAtVerts, !grepl('TxVP', species)) # remove fossils
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

T8 <- dplyr::select(PoAtVerts_wofossil_noNA, contains("Cen"), contains("c", ignore.case = FALSE), contains("specimen"), contains("species"))
T8 <- T8 %>% dplyr::rename(M1c = 3, M2c=4, M3c=5,M4c=6,M5c=7,M6c=8, M7c=9)

T12<- dplyr::select(PoAtVerts_wofossil_noNA, contains("d", ignore.case = FALSE), contains("specimen"), contains("species"))
T12 <- T12 %>% dplyr::rename(M1d = 1, M2d=2, M3d=3,M4d=4,M5d=5,M6d=6, M7d=7)

Sc <- dplyr::select(PoAtVerts_wofossil_noNA, contains("M", ignore.case = FALSE), contains("specimen"), contains("species"))



## PCA's ##

T1.pca <- prcomp(T1[c(1:7)], center = TRUE, scale = FALSE) # PCA
T4.pca <- prcomp(T4[c(1:7)], center = TRUE, scale = FALSE) # PCA
T8.pca <- prcomp(T8[c(1:9)], center = TRUE, scale = TRUE) # PCA
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
T8.rf <- randomForest(species ~ Cen_to_NeuAr+Cen_to_PoZy+ M1c + M2c+ M3c+M4c+M5c+M6c+M7c, data=T8_sub, importance=TRUE,proximity=TRUE)
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


### K Nearest neighbor ###:Non-parametric

library(caret)

#make KNN model using LOOCV to find optimal k

set.seed(123)

KNNmodel <- train(
  species ~., data = T8_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 6)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(T8_sub[,1:9]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == T8_sub$species) #overall accuracy

# assess accuracy per species
accKNN <- table(T8_sub$species,predicted.classes)
accKNN
diag(prop.table(accKNN, 1))

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


