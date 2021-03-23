## Load in data ##
library(readxl)
download.file(   "https://github.com/TIMAVID/Ambystoma/blob/master/Linear_data/Data/Amb_linear_data.xlsx?raw=true",   "Amb_linear_data.xlsx")
Amb_linear_data <- read_excel("Amb_linear_data.xlsx")


### ATLAS ###

# Tidy data #
library(dplyr)

Atlas <-  Amb_linear_data[c(3, 9:14, 55:56)] #select only relevant atlas measurements

Atlas_wofossil <- dplyr::filter(Atlas, !grepl('TxVP', species)) # remove fossils
Atlas_wofossil <- Atlas_wofossil %>% dplyr::rename(M1 = 2, M2=3, M3=4,M4=5,M5=6,M6=7) # remame columns for easier manipulation

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
autoplot(Atlas.pca_2, data = Atlas_wofossil_noTub, colour = 'species', frame = TRUE, label = TRUE) 

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
Atlas_wofossil_noTub_sub %>% cor_test(,1:6)

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

#pairwise comparisons between group levels with corrections for multiple testing
pairwise.perm.manova(Atlas_wofossil_noTub_sub[,1:6],Atlas_wofossil_noTub_sub$species,nperm=200) #needs more permutation but takes a long time
#or using euclidean distances
AtlasPPM<-pairwise.perm.manova(dist(Atlas_wofossil_noTub_sub[,1:6],"euclidean"),Atlas_wofossil_noTub_sub$species,nperm=999, progress = FALSE)
AtlasPPM


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
Atlascva=CVA(Atlas_wofossil_noTub_sub[,1:6], groups=Atlas_wofossil_noTub_sub$species, rounds = 10000, cv = TRUE)

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

# Plot Mahalahobis distances as dendrogram #

dendroS=hclust(Atlascva$Dist$GroupdistMaha)
dendroS$labels=levels(Atlas_wofossil_noTub_sub$species)
par(mar=c(6.5,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="",
     ylab='Mahalahobis distance')

### Random Forest ###:Non-parametric

library(randomForest)

Atlas.rf <- randomForest(species ~ M1 + M2 + M3 + M4 + M5 + M6, data=Atlas_wofossil_noTub_sub, importance=TRUE,
                         proximity=TRUE)
print(Atlas.rf)
rf_acc <- Atlas.rf$confusion
rf_acc <- 1-rf_acc[,11] # percent correct classification
rf_acc
# Look at variable importance
round(importance(Atlas.rf), 2)
varImpPlot(Atlas.rf)

### K Nearest neighbor ###:Non-parametric

library(tidyverse)
library(caret)

Atlas_wofossil_noTub_sub <- column_to_rownames(Atlas_wofossil_noTub_sub, var = "specimen_num")

#make KNN model using LOOCV to find optimal k
KNNmodel <- train(
  species ~., data = Atlas_wofossil_noTub_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 6)

plot(KNNmodel) # plot accuracy vs k
KNNmodel$bestTune # optimal k

predicted.classes <- KNNmodel %>% predict(Atlas_wofossil_noTub_sub[,1:6]) # predict class based on KNN model
head(predicted.classes)
mean(predicted.classes == Atlas_wofossil_noTub_sub$species) #overall accuracy

# assess accuracy per species
accKNN <- table(Atlas_wofossil_noTub_sub$species,predicted.classes)
accKNN
diag(prop.table(accKNN, 1))
