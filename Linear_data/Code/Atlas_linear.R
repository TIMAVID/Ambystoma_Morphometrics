## Load in data ##
library(readxl)
download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/Linear_data/Data/Amb_linear_data.xlsx?raw=true",   "Amb_linear_data.xlsx")
Amb_linear_data <- read_excel("Amb_linear_data.xlsx")


### ATLAS ###

# Tidy data #
library(dplyr)

Atlas <-  Amb_linear_data[c(1, 3, 9:14, 55:56)] #select only relevant atlas measurements
Atlas <- Atlas %>% dplyr::rename(M1 = 3, M2=4, M3=5,M4=6,M5=7,M6=8)# remame columns for easier manipulation

Atlas_wofossil <- dplyr::filter(Atlas, !grepl('41229*', species)) # remove fossils
Atlas_wofossil <- Atlas_wofossil[,-1]

Atlas_wofossil_noNA <- na.omit(Atlas_wofossil) # remove rows with N/A's

Atlas_wofossil_noNA$species<-as.factor(Atlas_wofossil_noNA$species)
Atlas_wofossil_noNA$species <- factor(Atlas_wofossil_noNA$species, levels = 
                                           c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                             "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")) # Reorder species
levels(Atlas_wofossil_noNA$species)

## PCA's ##
library(factoextra)
library(psych)

Atlas.pca <- prcomp(Atlas_wofossil_noNA[c(1:7)], center = TRUE, scale = TRUE) # PCA

# Summary stats #
summary(Atlas.pca)
sd <- Atlas.pca$sdev
loadings <- Atlas.pca$rotation
rownames(loadings) <- colnames(Atlas[c(1:7)])
scores <- Atlas.pca$x

# Show variance explained by PC's #
var <- sd^2
varPercent <- var/sum(var) * 100
barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')
abline(h=1/ncol(Atlas[c(1:7)])*100, col='red')

fviz_eig(Atlas.pca, xlab = "Principal Components")

#Show loadings #
library(devtools)
library(ggbiplot)

loadings
sqrt(1/ncol(Atlas[c(1:7)])) # cutoff for 'important' loadings

dev.new(height=7, width=7)
biplot(scores[, 1:2], loadings[, 1:2], cex=0.9)

ggbiplot(Atlas.pca, ellipse=FALSE, groups=Atlas_wofossil_noNA$species)

# Plot #
scores <- as.data.frame(scores)
scores$species <- Atlas_wofossil_noNA$species # reattach species
library(ggplot2)
library(grid)
library(gridExtra)
library(ggalt)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

percentage <- paste(colnames(scores), "(", paste(as.character(round(varPercent)), "%", " )", sep="") )

library(randomcoloR)
n <- 14
palette <- distinctColorPalette(n)

p<-ggplot(scores,aes(x=PC1,y=PC2,color=`species` ))
p<-p+geom_point(size =5)+theme + xlab(percentage[1]) + ylab(percentage[2]) +scale_color_manual(values = palette)
# p + stat_ellipse()
p
#alternative plot
fviz_pca_ind(Atlas.pca)
fviz_pca_ind(Atlas.pca, label="none", habillage=Atlas_wofossil_noNA$species,
             addEllipses=TRUE, ellipse.level=0.95, palette = palette)
#alternative plot 2
library(ggfortify)
autoplot(Atlas.pca, data = Atlas_wofossil_noNA, colour = 'species', frame = TRUE, label = TRUE) 


## No tuberculum interglenoideum ventral extent measurement ##
Atlas_wofossil_noTub <- Atlas_wofossil[,-1] 
Atlas_wofossil_noTub <- na.omit(Atlas_wofossil_noTub) # remove rows with N/A's

Atlas_wofossil_noTub$species<-as.factor(Atlas_wofossil_noTub$species)
Atlas_wofossil_noTub$species <- factor(Atlas_wofossil_noTub$species, levels = 
                                        c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                          "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")) # Reorder species

Atlas.pca_2 <- prcomp(Atlas_wofossil_noTub[c(1:6)], center = TRUE, scale = FALSE) # PCA
PC_scores <- as.data.frame(Atlas.pca_2$x)

percentage <- round(Atlas.pca_2$sdev / sum(Atlas.pca_2$sdev) * 100, 2)
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

autoplot(Atlas.pca_2, data = Atlas_wofossil_noTub, colour = 'species', frame = TRUE, label = TRUE) 

# PCA with fossils #

Atlas_fossil <- dplyr::filter(Atlas, grepl('41229*', species)) # fossils
Atlas_fossil <- subset(Atlas_fossil, select=-c(specimen_num, Specimen))


Atlas_fossil_complete <- na.omit(Atlas_fossil) # remove rows with N/A's

Amb_fossil_PCA <- predict(Atlas.pca_2, Atlas_fossil_complete[,2:7])
Fossil_PC_scores <- as.data.frame(Amb_fossil_PCA)

PC_scores <- cbind(PC_scores, species= Atlas_wofossil_noTub$species)
Fossil_PC_scores <- cbind(Fossil_PC_scores, species= Atlas_fossil_complete$species)

All_PC_scores <- rbind(PC_scores, Fossil_PC_scores) # create a new dataframe with the original PC scores and the PC scores of your fossil
tail(All_PC_scores)
pointsToLabel <- as.character(Atlas_fossil_complete$species)


species <- c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
             "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum")



pcaplot <- ggplot(data = All_PC_scores, mapping = aes(x = PC1, y = PC2,
                                                      col = species, label = species)) # creates the initial plot with datapoints color-coded and unique symbols by each species
pcaplot <- pcaplot + geom_encircle(expand=0, size = 2, data = All_PC_scores[!All_PC_scores$species %in% pointsToLabel,])+ theme_classic()
pcaplot <- pcaplot + 
  geom_text(aes(PC1, PC2, label = species), nudge_y = .003,
            check_overlap = FALSE, data = All_PC_scores[All_PC_scores$species %in% pointsToLabel,])+ geom_point(data = All_PC_scores)+ xlab(percentage[1]) + ylab(percentage[2])

pcaplot <- pcaplot + scale_color_manual(name = "Species", breaks = c(species),
                                        values=c("black", "black", "#D5E25E", "#AA47E3" ,"#8E7BD9" ,"#D2A6D5" ,"#7AA9D2" ,"#78DDD0", "#CAE1AE", "#D7A79D", "#DAB059", "#75E555", "#79E194",
                                                 "#CDDADD", "#DC5956", "#E363BB")) 
pcaplot


# Tuberculum interglenoideum plot

library(EnvStats)
Tub_dat <- Atlas_wofossil_noNA[c(1,9)]
  
ventral_extension_p <- ggplot(data = Tub_dat, aes(x = species, y = (tub_interglen_extension)))
ventral_extension_p <- ventral_extension_p + geom_boxplot(na.rm = TRUE)
ventral_extension_p <- ventral_extension_p + theme(axis.text.x = element_text(angle = 90))
ventral_extension_p <- ventral_extension_p + ylab("ventral extension (mm)") + stat_n_text() +theme_classic()
ventral_extension_p


### Statistical Tests ###

#Removing 'A.laterale|A.talpoideum|A.subsalsum|A.ordinarium' due to low sample sizes#

Atlas_wofossil_noTub_sub <- dplyr::filter(Atlas_wofossil_noTub, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))
Atlas_wofossil_noTub_sub$species <- factor(Atlas_wofossil_noTub_sub$species, levels = 
                                         c("A.gracile", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum",
                                           "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium")) # Reorder species

## Various ckecks for MANOVA ##

library(tidyverse)
library(ggpubr)
library(rstatix)
library(car)
library(broom)

# Atlas_wofossil_noTub_sub <- Atlas_wofossil_noTub_sub %>%
#   add_column(id = rownames(Atlas_wofossil_noTub_sub), .after = 8)

#Check sample sizes:PASS
Atlas_wofossil_noTub_sub %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

# Identify univariate outliers for each variable:FAIL
Atlas_wofossil_noTub_sub %>%
  dplyr::group_by(species) %>%
  identify_outliers(5) #input variable column here
Atlas_wofossil_noTub_sub %>%
  dplyr::group_by(species) %>%
  identify_outliers(6) #input variable column here
#...

#Detect multivariate outliers:PASS

mahalanobis_distance(Atlas_wofossil_noTub_sub[,1:6])

#Check univariate normality assumption:FAIL
Atlas_wofossil_noTub_sub %>%
  group_by(species) %>%
  shapiro_test(M1, M2,M3,M4,M5,M6) %>%
  arrange(variable)

#Check Multivariate normality:FAIL
Atlas_wofossil_noTub_sub %>%
  dplyr::select(,1:6) %>%
  mshapiro_test()

#Identify multicollinearity:FAIL
Atlas_wofossil_noTub_sub %>% rstatix::cor_test(,1:6)

# PROBLEM!!:Absence of multicollinearity. The dependent (outcome) variables cannot be too correlated to each other. No correlation should be above r = 0.90 [Tabachnick and Fidell (2012)}.
cor(Atlas_wofossil_noTub_sub[,1:6])

#MANOVA# :*Failed multiple checks
Atlas.man <- manova(cbind(M1,M2,M3,M4,M5,M6) ~ species, data = Atlas_wofossil_noTub_sub)
summary(Atlas.man)
summary.aov(Atlas.man)

## Permutation MANOVA ## *Overcome failed checks
library(RVAideMemoire)
require(vegan)
#Permutational Multivariate Analysis of Variance Using Distance Matrices
adonis(Atlas_wofossil_noTub_sub[,1:6]~species,data=Atlas_wofossil_noTub_sub,method="euclidean") 
# ?adonis

#pairwise comparisons between group levels with corrections for multiple testing
pairwise.perm.manova(Atlas_wofossil_noTub_sub[,1:6],Atlas_wofossil_noTub_sub$species,nperm=50) #needs more permutation but takes a long time
#or using euclidean distances
AtlasPPM<-pairwise.perm.manova(dist(Atlas_wofossil_noTub_sub[,1:6],"euclidean"),Atlas_wofossil_noTub_sub$species,nperm=999, progress = FALSE)
AtlasPPM
t <- AtlasPPM$p.value
t <-round(t, digits = 3)
write.table(t, file = "Atlas linear PW", sep = ",", quote = FALSE, row.names = T)

# tuberculum interglenoideum measurement only
Atlas_wofossil_Tub_only <- as.data.frame(Atlas_wofossil[c(1,8:9)]) 
Atlas_wofossil_Tub_only <- na.omit(Atlas_wofossil_Tub_only) # remove rows with N/A's
Atlas_wofossil_Tub_only$species<-as.factor(Atlas_wofossil_Tub_only$species)
Atlas_wofossil_Tub_only$species <- factor(Atlas_wofossil_Tub_only$species, levels = 
                                             c("A.gracile", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum",
                                               "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium")) # Reorder species


Atlas_wofossil_Tub_only <- dplyr::filter(Atlas_wofossil_Tub_only, !grepl('A.laterale|A.talpoideum|A.subsalsum|A.ordinarium', species))


#Permutational Anova
set.seed(123)
perm.anova(Atlas_wofossil_Tub_only$tub_interglen_extension ~ Atlas_wofossil_Tub_only$species, nperm=1000)
?perm.anova
#pairwise comparisons between group levels with corrections for multiple testing
library(rcompanion)

PT <- pairwisePermutationTest(tub_interglen_extension ~   species,
                             data   = Atlas_wofossil_Tub_only,
                             method = "fdr")
PT

t<- pairwise.perm.t.test(Atlas_wofossil_Tub_only$tub_interglen_extension,Atlas_wofossil_Tub_only$species,nperm=999,progress = FALSE)
t <- t$p.value
t <- round(t, 3)
write.table(t, file = "Tub PW", sep = ",", quote = FALSE, row.names = T)
### DFA ###

#PROBLEM:Check Multivariate normality:FAIL
mqqnorm(Atlas_wofossil_noTub_sub[,1:6], main = "Multi-normal Q-Q Plot")

#DFA# With MASS

library(MASS)
AtlasLDA <- lda(species ~ M1 + M2 + M3 + M4 + M5 + M6, data=Atlas_wofossil_noTub_sub, CV = FALSE) #DFA no jacknife
AtlasLDA

AtlasLDA_jack <- lda(species ~ M1 + M2 + M3 + M4 + M5 + M6, data=Atlas_wofossil_noTub_sub, CV = TRUE) #DFA with jacknife
AtlasLDA_jack

# Assess the accuracy of jacknife #

accAtlasLDA <- table(Atlas_wofossil_noTub_sub$species, AtlasLDA_jack$class)
accAtlasLDA
diag(prop.table(accAtlasLDA, 1))
sum(accAtlasLDA[row(accAtlasLDA) == col(accAtlasLDA)]) / sum(accAtlasLDA)

#DFA# With MORPHO

library(Morpho)
Atlascva <- CVA(Atlas_wofossil_noTub_sub[,1:6], groups=Atlas_wofossil_noTub_sub$species, rounds = 0, cv = TRUE)

barplot(Atlascva$Var[,2]) # Variance explained by the canonical roots

# get the typicality probabilities and resulting classifications
# all specimens with a probability of < 0.01 as outliers (assigned to no class)
typprobs <- typprobClass(Atlascva$CVscores,groups=Atlas_wofossil_noTub_sub$species, outlier = 0.01, cv = TRUE)
print(typprobs)

# Assess the accuracy of jacknife #

accJack <- table(Atlascva$groups, Atlascva$class)
accJack
diag(prop.table(accJack, 1))
sum(accJack[row(accJack) == col(accJack)]) / sum(accJack)

# Plot first two DF axes #

AT_cva <- data.frame(Atlascva$CVscores, species = Atlascva$groups)

ggplot(AT_cva, aes(CV.1, CV.2)) +
  geom_point(size =5,aes(color = species)) + theme_classic() + scale_color_brewer(palette="Paired")
#alternative plot
plot(Atlascva$CVscores, col=Atlas_wofossil_noTub_sub$species, pch=as.numeric(Atlas_wofossil_noTub_sub$species), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(Atlascva$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(Atlascva$Var[2,2],1),"%")))
text(Atlascva$CVscores, as.character(Atlas_wofossil_noTub_sub$species), col=as.numeric(Atlas_wofossil_noTub_sub$species), cex=.7)

# DFA Fossil classification #

fossil_CVA_scores <- predict(Atlascva, as.matrix(Atlas_fossil_complete[,2:7]))

fossil_class <- classify(Atlascva, cv = FALSE, newdata = as.matrix(Atlas_fossil_complete[,2:7]))
fossil_class$class
fossil_class$posterior

plot(Atlascva$CVscores, col=Atlas_wofossil_noTub_sub$species, pch=as.numeric(Atlas_wofossil_noTub_sub$species), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(Atlascva$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(Atlascva$Var[2,2],1),"%")))

text(fossil_CVA_scores, as.character(Atlas_fossil_complete$species), cex=.7)
text(Atlascva$CVscores, as.character(Atlas_wofossil_noTub_sub$species), col=as.numeric(Atlas_wofossil_noTub_sub$species), cex=.7)

# Plot Mahalahobis distances as dendrogram #

dendroS=hclust(Atlascva$Dist$GroupdistMaha)
dendroS$labels=levels(Atlas_wofossil_noTub_sub$species)
par(mar=c(6.5,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="",
     ylab='Mahalahobis distance')

### Random Forest ###:Non-parametric

library(randomForest)
set.seed(123)
Atlas.rf <- randomForest(species ~ M1 + M2 + M3 + M4 + M5 + M6, data=Atlas_wofossil_noTub_sub, importance=TRUE,
                         proximity=TRUE)
print(Atlas.rf)
rf_acc <- Atlas.rf$confusion
rf_acc <- 1-rf_acc[,11] # percent correct classification
rf_acc

t <- rf_acc
t <-round(t, digits = 2)
write.table(t, file = "Atlas linear KNNAC", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf$predicted == Atlas_wofossil_noTub_sub$species) #overall accuracy

# Look at variable importance
round(importance(Atlas.rf), 2)
varImpPlot(Atlas.rf)

# with variable M4 removed
set.seed(123)
Atlas.rf_M1 <- randomForest(species ~ M1 + M2 + M3 + M5 + M6, data=Atlas_wofossil_noTub_sub, importance=TRUE,
                         proximity=TRUE, replace=FALSE)
print(Atlas.rf_M1)
rf_acc_M1 <- Atlas.rf_M1$confusion
rf_acc_M1 <- 1-rf_acc_M1[,11] # percent correct classification
rf_acc_M1
mean(Atlas.rf_M1$predicted == Atlas_wofossil_noTub_sub$species) #overall accuracy


# Predict fossils

y_pred = predict(Atlas.rf, newdata = Atlas_fossil_complete[,2:7])
y_pred

t <- y_pred
t <- as.data.frame(t)
t <-round(t, digits = 2)
write.table(t, file = "Atlas linear RFAC", sep = ",", quote = FALSE, row.names = T)

### K Nearest neighbor ###:Non-parametric

library(tidyverse)
library(caret)

Atlas_wofossil_noTub_sub <- column_to_rownames(Atlas_wofossil_noTub_sub, var = "specimen_num")

#make KNN model using LOOCV to find optimal k

set.seed(123)

KNNmodel <- train(
  species ~., data = Atlas_wofossil_noTub_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 10)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(Atlas_wofossil_noTub_sub[,1:6]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == Atlas_wofossil_noTub_sub$species) #overall accuracy

# assess accuracy per species
accKNN <- table(Atlas_wofossil_noTub_sub$species,predicted.classes)
accKNN
t <- diag(prop.table(accKNN, 1))
t <-round(t, digits = 2)
write.table(t, file = "Atlas linear KNNAC", sep = ",", quote = FALSE, row.names = T)

# Fossil predictions #

library(class)
KnnTestPrediction_k9 <- knn(Atlas_wofossil_noTub_sub[,1:6], Atlas_fossil_complete[,2:7],
                            Atlas_wofossil_noTub_sub$species, k=9, prob=TRUE)
KnnTestPrediction_k9

KnnTestPrediction_k5 <- knn(Atlas_wofossil_noTub_sub[,1:6], Atlas_fossil_complete[,2:7],
                            Atlas_wofossil_noTub_sub$species, k=5, prob=TRUE)
KnnTestPrediction_k5

KnnTestPrediction_k3 <- knn(Atlas_wofossil_noTub_sub[,1:6], Atlas_fossil_complete[,2:7],
                            Atlas_wofossil_noTub_sub$species, k=3, prob=TRUE)
KnnTestPrediction_k3


### Model selection (multinominal regression) 

library(glmulti)

library(nnet)

multinom.glmulti <- function(formula, data, ...)
  multinom(formula, data, ...)

res <- glmulti(species ~ .,level=1, data=Atlas_wofossil_noTub_sub, report = FALSE, plotty = FALSE,fitfunction=multinom.glmulti, method = "h", crit="aic", confsetsize=16)

print(res)

plot(res)

top <- weightable(res)
top



multinom_model1 <- nnet::multinom(species ~ M1 + M3 + M4 + M5 + M6, data = Atlas_wofossil_noTub_sub, maxit=1000)
multinom_model2 <- nnet::multinom(species ~ M1 + M3 + M4 + M5, data = Atlas_wofossil_noTub_sub, maxit=1000)
multinom_model3 <- nnet::multinom(species ~ M1 + M2 + M3 + M4 + M5 + M6, data = Atlas_wofossil_noTub_sub, maxit=1000)
multinom_model4 <- nnet::multinom(species ~ M1 + M2 + M3 + M5 + M6, data = Atlas_wofossil_noTub_sub, maxit=1000)
multinom_model5 <- nnet::multinom(species ~ M1 + M2 + M3 + M6, data = Atlas_wofossil_noTub_sub, maxit=1000)
multinom_model6 <- nnet::multinom(species ~ M2 + M3, data = Atlas_wofossil_noTub_sub, maxit=1000)


modelacc <- function(model_name,test.data, species){
  predicted.classes <- model_name %>% predict(test.data)
  print(mean(predicted.classes == species))
}

modelacc(multinom_model6, Atlas_wofossil_noTub_sub, Atlas_wofossil_noTub_sub$species)


predicted.classes <- multinom_model6 %>% predict(Atlas_fossil_complete[,2:7])









# KNN with top models

set.seed(123)

KNNmodel_1 <- train(
  species ~M1 + M3 + M4 + M5 + M6, data = Atlas_wofossil_noTub_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 10)

plot(KNNmodel_1) # plot accuracy vs k
KNNmodel_1$bestTune # optimal k

predicted.classes_M1 <- KNNmodel_1 %>% predict(Atlas_wofossil_noTub_sub[c(1,3:6)]) # predict class based on KNN model
head(predicted.classes_M1)
mean(predicted.classes_M1 == Atlas_wofossil_noTub_sub$species) #overall accuracy


set.seed(123)

KNNmodel_2 <- train(
  species ~M1 + M3 + M4 + M5, data = Atlas_wofossil_noTub_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 10)

plot(KNNmodel_2) # plot accuracy vs k
KNNmodel_2$bestTune # optimal k

predicted.classes_M2 <- KNNmodel_2 %>% predict(Atlas_wofossil_noTub_sub[c(1,3:5)]) # predict class based on KNN model
head(predicted.classes_M2)
mean(predicted.classes_M2 == Atlas_wofossil_noTub_sub$species) #overall accuracy


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
Amb_species<-unique(Atlas_wofossil_noTub_sub$species)
tips<-tree$tip.label
ii<-sapply(Amb_species,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
plotTree(tree,ftype="i")
#Tree did not include A.mavortium so I lumped that species with A.tigrinum
Atlas_wofossil_noTub_sub$species<-gsub("A.mavortium", "A.tigrinum", Atlas_wofossil_noTub_sub$species, fixed = TRUE)
Atlas_wofossil_noTub_sub$species<-as.factor(Atlas_wofossil_noTub_sub$species)
Atlas_wofossil_noTub_sub$species <- factor(Atlas_wofossil_noTub_sub$species, levels = 
                                             c("A.gracile", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum",
                                               "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium")) # Reorder species


#Preformed a group PCA
library(Morpho)
library(geomorph)
gpca <- groupPCA(Atlas_wofossil_noTub_sub[,1:6], Atlas_wofossil_noTub_sub$species, rounds=0)
plot(gpca$groupmeans)
?groupPCA
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




