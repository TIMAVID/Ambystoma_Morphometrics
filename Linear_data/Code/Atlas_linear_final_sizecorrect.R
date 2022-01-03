## Load in data ##
# Read in csv file of data
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/Linear_data/Data/Amb_linear_data_final.csv")
Amb_linear_data <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(Amb_linear_data)

### ATLAS ###

# Tidy data #
library(dplyr)

Atlas <-  Amb_linear_data[c(2,3, 8:13,50, 54:55)] #select only relevant atlas measurements

# tub_interglen_extension = ventral extent of the atlantal cotyles below the tuberculum interglenoideum
# 1 = mid-ventral length of the atlas
# 2 = width between the atlantal cotyles
# 3 = height of the atlantal cotyle
# 4 = width between the postzygapophyses
# 5 = height from the mid-posteroventral margin of the centrum to the dorsal apex of the neural arch
# 6 = width at the posterior end of the centrum






# REMOVE FOSSILS AND TUBERCULUM INTERGLENOIDEUM MEASUREMENT #---------------------------------
Atlas_wofossil <- dplyr::filter(Atlas, !grepl('41229*', species)) # remove fossils
row.names(Atlas_wofossil) <- Atlas_wofossil$specimen_num
Atlas_wofossil$species <- factor(Atlas_wofossil$species, levels = 
                                   c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                     "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels

Atlas_wofossil_noTub <- subset(Atlas_wofossil, select=-c(tub_interglen_extension, Cotyle_height)) # no tuberculum interglenoideum extension measurement
Atlas_wofossil_noTub <- na.omit(Atlas_wofossil_noTub) # remove rows with N/A's
Atlas_wofossil_noTub <- subset(Atlas_wofossil_noTub, select=-c(specimen_num))
species <- Atlas_wofossil_noTub$species


## SUBSET DATA FOR FOSSILS ONLY ##---------------------------------

Atlas_fossil <- dplyr::filter(Atlas, grepl('41229*', species)) # fossils
Atlas_fossil <- subset(Atlas_fossil, select=-c(SVL_P, specimen_num)) # no tuberculum interglenoideum extension measurement
Atlas_fossil <- na.omit(Atlas_fossil) # remove rows with N/A's
row.names(Atlas_fossil) <- Atlas_fossil$species
Atlas_fossil_noTub <- subset(Atlas_fossil, select=-c(tub_interglen_extension, Cotyle_height))
fossils <- Atlas_fossil_noTub$species


# EFFECT OF SIZE ##---------------------------------

Gms<- apply(Atlas_wofossil_noTub[,1:6], 1, function(x) exp(mean(log(x))))

total <- data.frame(SVL_p=Atlas_wofossil_noTub$SVL_P, Gms, species=Atlas_wofossil_noTub$species)

fit1 <- lm(Gms ~ log(SVL_p), data = total)
summary(fit1)

speciescolors <- c("#666600", "#C9D42D" ,"#42CEBC" ,"#F2AB1F" ,"#864ED0" ,"#261ACE", "#086AFD", "#08FD6A", "#0C8B3F", "#E50CF5", "#FF5E00","#FF0000", "#FF6A6A", "#D5930F", "#9E1616", "#000000", "#000000" )

ggplot(total, aes(log(SVL_p), Gms, color=species)) + 
  geom_point() +  stat_smooth(method = "lm", col = "red") +
  theme_classic()+ coord_fixed(ratio = 1)+
  scale_color_manual(name = "Species", breaks=levels(total$species), values=c(speciescolors))


# FUNCTIONS TO CREATE LINEAR MODELS FOR EXAMINING ALLOMETRY------------
plot.alom.species<-function(variables,data){
  require(ggplot2)
  require(gridExtra)
  figs<-lapply(variables, function(x) {
    ggplot(data = data,
           aes(log(SVL_P), log(get(x)))) + geom_point(aes(color=species))+
      scale_colour_manual(values=speciescolors)+ggtitle(x)+
      theme_classic(base_size = 8)+ ylab("log(Measurement)")+
      geom_smooth(method='lm')+theme(legend.position="none")+ coord_fixed(ratio = 1)+
      geom_abline(intercept=seq(-100, 100, 1),
                  slope=1,
                  colour="red",linetype='dashed')
  })
  do.call(grid.arrange, c(figs, ncol=3, top = deparse(substitute(data))))
  
  models <- lapply(variables, function(x) { #function to perform linear regression on all measurements for each dataset
    lm(substitute(log(i) ~ log(SVL_P), list(i = as.name(x))), data = data)
  })
  names(models) <- variables
  (sum <- (lapply(models, summary)))
  b<-(lapply(sum, function (x)x$coefficients[1]))
  m<-(lapply(sum, function (x)x$coefficients[2]))
  R<-(lapply(sum, function (x)x$adj.r.squared))
  P<-(lapply(sum, function (x)x$coefficients[2,4])) #may need to p.adjust
  out<-list(m,b,R,P,models)
  names(out)<-c("slope","Y-intercept","adj.R-squared","P-value","models")
  out <- do.call(rbind, out)
  out <- t(out)
  out <- as.data.frame(out)
  return(out)
}
varlist <- names(Atlas_wofossil_noTub)[1:6] #all the different measurements

Amby_atlas_alom <-plot.alom.species(varlist,Atlas_wofossil_noTub[,1:7]) #liner models of all measurements


# REMOVE EFFECT OF "SIZE" ##---------------------------------

library(tidyverse)

geometric_mean <- function(x){
  exp(sum(log(x), na.rm = TRUE) / length(x))
}

Shape.varb.calc <- function(x){
  log(x / geometric_mean(x))
}

Atlas_wofossil_noTub <- as.data.frame(t(apply(Atlas_wofossil_noTub[,1:6], 1, Shape.varb.calc)))
Atlas_wofossil_noTub$species <- species

Atlas_fossil_noTub <- as.data.frame(t(apply(Atlas_fossil_noTub[,1:6], 1, Shape.varb.calc)))
Atlas_fossil_noTub$species <- fossils


ggplot(total, aes(Gms, Atlas_wofossil_noTub[,2], color=species)) + 
  geom_point() +  stat_smooth(method = "lm", col = "red") + theme_bw()+
  scale_color_manual(name = "Species", breaks=levels(total$species), values=c(speciescolors))


## PRINCIPAL COMPONENT ANALYSIS ##---------------------------------

Atlas.pca <- prcomp(Atlas_wofossil_noTub[c(1:6)], center = TRUE, scale = FALSE) # PCA
PC_scores <- as.data.frame(Atlas.pca$x)
PC_scores <- cbind(PC_scores, species= species)
atlasloadings <- data.frame(Variables = rownames(Atlas.pca$rotation), Atlas.pca$rotation)# Extract loadings of the variables

percentage <- round(Atlas.pca$sdev^2 / sum(Atlas.pca$sdev^2) * 100, 2)# find percentage variance explained by PC's
percentage <- paste( colnames(PC_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

# PROJECT FOSSIL DATA #
Amb_fossil_PCA <- predict(Atlas.pca, Atlas_fossil_noTub[,1:6])
Fossil_PC_scores <- as.data.frame(Amb_fossil_PCA)
Fossil_PC_scores <- cbind(Fossil_PC_scores, species= fossils)
Fossil_PC_scores <- droplevels(Fossil_PC_scores)

All_PC_scores <- rbind(PC_scores, Fossil_PC_scores) # create a new dataframe with the original PC scores and the PC scores of your fossil
tail(All_PC_scores)

# PLOT #
speciescolors <- c("#666600", "#C9D42D" ,"#42CEBC" ,"#F2AB1F" ,"#864ED0" ,"#261ACE", "#086AFD", "#08FD6A", "#0C8B3F", "#E50CF5", "#FF5E00","#FF0000", "#FF6A6A", "#D5930F", "#9E1616", "black", "black" )
speciesshapes <- c(rep(16,15), rep(18,30))

library(ggplot2)
library(ggforce)
p<-ggplot(All_PC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_segment(data = atlasloadings, aes(x = 0, y = 0, xend = (.83*PC1),
                                      yend = (.8*PC2)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +annotate("text", x = (atlasloadings$PC1), y = (atlasloadings$PC2),
                                          label = atlasloadings$Variables) +
  geom_point(size =3)+ xlab(percentage[1]) + ylab(percentage[2]) + coord_fixed()+
  scale_shape_manual(values = c(speciesshapes), guide = 'none') + theme_classic() + theme(legend.position = "none")+ 
  scale_color_manual(name = "Species", breaks=levels(All_PC_scores$species), values=c(speciescolors))
p

library(factoextra)
fviz_pca_var(Atlas.pca,
             #col.var = "contrib", # Color by contributions to the PC
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# 3D plot
library(plotly)
plot_ly(x=All_PC_scores$PC1, y=All_PC_scores$PC2, z=All_PC_scores$PC3, type="scatter3d", mode="markers", color=All_PC_scores$species)





# *removed A. subsalsum and A. ordinarium* see code above for removal process
Atlas_wofossil_noTub_sub <- dplyr::filter(Atlas_wofossil_noTub, !grepl('A.subsalsum|A.ordinarium', species))
Atlas_wofossil_noTub_sub$species <- factor(Atlas_wofossil_noTub_sub$species, levels = 
                                             c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                               "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.velasci")) # Reorder species levels

# Sample sizes #
Atlas_wofossil_noTub_sub %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())







### K NEAREST NEIGHBOR   ###:Non-parametric---------------------------------
library(caret)


# LOOCV WITH REPLICATION
library(foreach)
library(doParallel)
ncore <- detectCores()
registerDoParallel(cores=ncore)

set.seed(123)

runs <- 100

system.time({
  fish <- foreach(icount(runs)) %dopar% {
    train(species~ .,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = Atlas_wofossil_noTub_sub)$results
  }
})

fish <- map_dfr(fish,`[`, c("k", "Accuracy", "Kappa"))
kspecies <- fish %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kspecies <- kspecies[1,]



library(class)
set.seed(123)
predicted.classes <- train(species~ .,
                           method     = "knn",
                           tuneGrid   = expand.grid(k = kspecies),
                           trControl  = trainControl(method  = "LOOCV"),
                           metric     = "Accuracy",
                           data       = Atlas_wofossil_noTub_sub)$pred # predict class based on KNN model
mean(predicted.classes$pred == predicted.classes$obs) #overall accuracy

accKNN <- table(predicted.classes$obs,predicted.classes$pred)
accKNN
# t <- diag(prop.table(accKNN, 1))
# t <-round(t, digits = 2)
# write.table(t, file = "Atlas linear species KNNAC size corrected", sep = ",", quote = FALSE, row.names = T)

# FOSSIL CLASSIFICATION #

KnnTestPrediction_k6 <- knn(Atlas_wofossil_noTub_sub[,1:6], Atlas_fossil_noTub[,1:6],
                            Atlas_wofossil_noTub_sub$species, k=6, prob=TRUE)
KnnTestPrediction_k6
# write.table(KnnTestPrediction_k6, file = "Atlas fossil KNN size corrected", sep = ",", quote = FALSE, row.names = T)





### RANDOM FOREST CLASSIFICATION ###:Non-parametric---------------------------------
library(randomForest)
set.seed(123)
Atlas.rf <- randomForest(species ~ ., data=Atlas_wofossil_noTub_sub, importance=TRUE)
print(Atlas.rf)
rf_acc <- Atlas.rf$confusion
rf_acc <- 1-rf_acc[,14] # percent correct classification
rf_acc
t <- rf_acc
t <-round(t, digits = 2)
write.table(t, file = "size corrected Atlas linear RFAC", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf$predicted == Atlas_wofossil_noTub_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred = predict(Atlas.rf, newdata = Atlas_fossil_noTub[,1:6])
y_pred
# write.table(y_pred, file = "size corrected Atlas fossil RF", sep = ",", quote = FALSE, row.names = T)






# AMBYSTOMA CLADE CLASSIFICATION #---------------------------------

species <- Atlas_wofossil_noTub_sub$species
clades <- dplyr::recode(species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
Atlas_wofossil_noTub_sub_clade <- data.frame(Atlas_wofossil_noTub_sub[,1:6],clades=clades)

#KNN#
set.seed(123)

system.time({
  fish2 <- foreach(icount(runs)) %dopar% {
    train(clades~ .,
          method     = "knn",
          tuneGrid   = expand.grid(k = 1:17),
          trControl  = trainControl(method  = "LOOCV"),
          metric     = "Accuracy",
          data       = Atlas_wofossil_noTub_sub_clade)$results
  }
})


fish2 <- map_dfr(fish2,`[`, c("k", "Accuracy", "Kappa"))
kclade <- fish2 %>% 
  filter(Accuracy == max(Accuracy)) %>% # filter the data.frame to keep row where Accuracy is maximum
  select(k) # select column k
kclade <- kclade[1,]

set.seed(123)
predicted.clades <- train(clades~ .,
                          method     = "knn",
                          tuneGrid   = expand.grid(k = kclade),
                          trControl  = trainControl(method  = "LOOCV"),
                          metric     = "Accuracy",
                          data       = Atlas_wofossil_noTub_sub_clade)$pred # predict class based on KNN model
mean(predicted.clades$pred == predicted.clades$obs) #overall accuracy

accKNN_clades <- table(predicted.clades$obs,predicted.clades$pred)
accKNN_clades
t <- diag(prop.table(accKNN_clades, 1))
t <-round(t, digits = 2)
write.table(t, file = "size corrected Atlas linear clades KNNAC", sep = ",", quote = FALSE, row.names = T)


# FOSSIL CLADE CLASSIFICATION #

KnnTestPrediction_k9 <- knn(Atlas_wofossil_noTub_sub_clade[,1:6], Atlas_fossil_noTub[,1:6],
                            Atlas_wofossil_noTub_sub_clade$clades, k=9, prob=TRUE)
KnnTestPrediction_k9
write.table(KnnTestPrediction_k9, file = "size corrected Atlas fossil clades KNN", sep = ",", quote = FALSE, row.names = T)






# RANDOM FOREST WITH CLADES #---------------------------------
set.seed(123)
Atlas.rf_clades <- randomForest(clades ~ ., data=Atlas_wofossil_noTub_sub_clade, importance=TRUE)
print(Atlas.rf_clades)
rf_acc_clades <- Atlas.rf_clades$confusion
rf_acc_clades <- 1-rf_acc_clades[,9] # percent correct classification
rf_acc_clades
t <- rf_acc_clades
t <-round(t, digits = 2)
write.table(t, file = "size corrected Atlas linear RFACclades", sep = ",", quote = FALSE, row.names = T)

mean(Atlas.rf_clades$predicted == Atlas_wofossil_noTub_sub_clade$clades) #overall accuracy

# FOSSIL CLADE CLASSIFICATION #
y_pred_clade = predict(Atlas.rf_clades, newdata = Atlas_fossil_noTub[,1:6])
y_pred_clade
write.table(y_pred_clade, file = "size corrected Atlas fossil RFclades", sep = ",", quote = FALSE, row.names = T)


# MEASUREMENT RELATIVE IMPORTANCE BASED ON RF #---------------------------------

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
          axis.text.y  = element_text(size = 15, color = "black")) + theme_classic()
  return(p)
}

create_crfplot(rf3, conditional = TRUE)


