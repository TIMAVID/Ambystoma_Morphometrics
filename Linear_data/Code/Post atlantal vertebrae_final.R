## Load in data ##
# Read in csv file of data
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
PoAtVerts_wofossil <- subset(PoAtVerts_wofossil, select=-c(specimen_num))
PoAtVerts_wofossil <- droplevels(PoAtVerts_wofossil)
PoAtVerts_wofossil$species <- factor(PoAtVerts_wofossil$species, levels = 
                                            c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                              "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.mavortium", "A.ordinarium", "A.subsalsum", "A.velasci")) # Reorder species levels




# SUBSET DATA FOR EACH VERTEBRA #
T1 <- dplyr::select(PoAtVerts_wofossil, contains("a", ignore.case = FALSE),  contains("species")) %>% drop_na()

T4 <- dplyr::select(PoAtVerts_wofossil, contains("b", ignore.case = FALSE), contains("species")) %>% drop_na()

T8 <- dplyr::select(PoAtVerts_wofossil, contains("c", ignore.case = FALSE), contains("species")) %>% drop_na()

T8_extension <- dplyr::select(PoAtVerts_wofossil, contains("Cen"), contains("species")) %>% drop_na()

T12<- dplyr::select(PoAtVerts_wofossil, contains("d", ignore.case = FALSE), contains("species")) %>% drop_na()

Sc <- dplyr::select(PoAtVerts_wofossil, num_range("X", 14:20), contains("species")) %>% drop_na()
Sc.pca <- prcomp(Sc[c(1:7)], center = TRUE, scale = FALSE) # PCA






### CALCULATE ZYGAPOPHYSEAL RATIOS ###
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

speciescolors <- c("#666600", "#999900" ,"#99FFFF" ,"#FF9933" ,"#B266FF" ,"#1139BC", "#003BFD", "#00FD33", "#26B543", "#E50CF5", "#FF0909","#D71D1D", "#EC6A6A", "#8B0B0B", "#A54E4E")

ZyPlot <- ggplot(ZySummary, aes(x=variable, y=mean, group=species, color=species)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.1, position=position_dodge(0.15)) +
  geom_line() + geom_point()+
  scale_color_manual(name = "Species", values=c(speciescolors)) + theme_minimal() + xlab("Vertebra") + ylab("Zygapoheseal ratio") 
ZyPlot <- ZyPlot + facet_wrap(~species, ncol = 15)  # wrap data 'by' family into 4 columns
ZyPlot + scale_x_discrete(breaks=c("Zyratioa","Zyratiob","Zyratioc", "Zyratiod", "ZyratioSc"),
                           labels=c("T1", "T4", "T8", "T12", "SC"))
ZyPlot <- ZyPlot + theme(legend.position = "none")






### CALCULATE CENTRUM RATIOS ###

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
                          labels=c("T1", "T4", "T8", "T12", "SC"))






## Phylogenetic signal ##

# Load in data #
require(phytools)

f3 <- curl("https://raw.githubusercontent.com/TIMAVID/Ambystoma/master/GMM/Data/Amb_species")

# Read in tree

tree <- read.newick(f3)  #tree from Williams et al. 2013

par(mar=c(1,1,1,1))
tree$tip.label<-gsub("^", "A.", tree$tip.label)
plot(tree)

#Subset tree to include only GMM species
Amb_species<-unique(T8$species)
tips<-tree$tip.label
ii<-sapply(Amb_species,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
plotTree(tree,ftype="i")

#Tree did not include A.mavortium so I lumped that species with A.tigrinum
T8_sub_tig<- ZySummary
T8_sub_tig$species<-gsub("A.mavortium", "A.tigrinum", T8_sub_tig$species, fixed = TRUE)
T8_sub_tig$species<-gsub("A.subsalsum", "A.ordinarium", T8_sub_tig$species, fixed = TRUE)
T8_sub_tig$species<-gsub("A.velasci", "A.ordinarium", T8_sub_tig$species, fixed = TRUE)
T8_sub_tig$species<-as.factor(T8_sub_tig$species)
T8_sub_tig$species <- factor(T8_sub_tig$species, levels = 
                               c("A.gracile","A.talpoideum", "A.maculatum", "A.macrodactylum","A.opacum","A.jeffersonianum","A.laterale",
                                 "A.mabeei","A.texanum","A.annulatum","A.tigrinum","A.ordinarium")) # Reorder species
T8_sub_tig <- T8_sub_tig %>% 
  dplyr::filter(grepl("Zyratioc", variable)) %>% 
  dplyr::group_by(species) %>%
  dplyr::summarise(Zymean = mean(mean))

library(geomorph)
test <- (T8_sub_tig$mean)
names <- T8_sub_tig$species
setNames(test, T8_sub_tig$species)

test$dimnames[[1]] <- as.character(names)

physignal(T8_sub_tig$mean, tree, print.progress = F, iter = 999)


names(T8_sub_tig$mean)













# MAKE REGIONAL VERTEBRAE GROUPING SUBDATASETS "ANTERIOR, MIDDLE, POSTERIOR" #

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





# MAKE FOSSIL VERTEBRAE GROUPINGS BASED ON ESTIMATED POSITION IN VERTEBRAL COLUMN #

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






## PRINCIPAL COMPONENT ANALYSES WITH FOSSILS ##

# Anterior trunk fossils

TrunkFossilAnt_PCA <- predict(TrunkAnt.pca, TrunkFossilAnt[,1:7])
TrunkFossilAnt_PC_scores <- as.data.frame(TrunkFossilAnt_PCA)

FossilAnt_PC_scores <-data.frame(TrunkFossilAnt_PCA, species = TrunkFossilAnt$species)

Antscores <-data.frame(TrunkAnt.pca$x, species = TrunkAnt$species)

All_AntPC_scores <- (rbind(Antscores, FossilAnt_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of the fossils

# PLOT #

fossilcolors <- grDevices::gray.colors(54, start = 0, end = 0)

library(ggplot2)
library(ggforce)

percentage_ant <- round(TrunkAnt.pca$sdev / sum(TrunkAnt.pca$sdev) * 100, 2)# find percentage variance explained by PC's
percentage_ant <- paste( colnames(Antscores), "(", paste( as.character(percentage_ant), "%", ")", sep="") )

Ant_plot<-ggplot(All_AntPC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_ant[1]) + ylab(percentage_ant[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkAnt$species), values=c(speciescolors, fossilcolors)) + theme_classic() + ggtitle("Anterior Vertebrae")
Ant_plot




# Middle trunk fossils

TrunkFossilMID_PCA <- predict(TrunkMid.pca, TrunkFossilT8[,1:7])
TrunkFossilMID_PC_scores <- as.data.frame(TrunkFossilMID_PCA)

MIDscores <-data.frame(TrunkMid.pca$x, species= (TrunkMid$species))

FossilMID_PC_scores <-data.frame(TrunkFossilMID_PCA, species= TrunkFossilT8$species)

All_MIDPC_scores <- (rbind(MIDscores, FossilMID_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil

percentage_mid <- round(TrunkMid.pca$sdev / sum(TrunkMid.pca$sdev) * 100, 2)# find percentage variance explained by PC's
percentage_mid <- paste( colnames(MIDscores), "(", paste( as.character(percentage_mid), "%", ")", sep="") )

# PLOT #
Mid_plot<-ggplot(All_MIDPC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_mid[1]) + ylab(percentage_mid[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkMid$species), values=c(speciescolors, fossilcolors)) + theme_classic() + ggtitle("Middle Vertebrae")
Mid_plot




# Posterior trunk fossils

TrunkFossilPOST_PCA <- predict(TrunkPost.pca, TrunkFossilT12[,1:7])
TrunkFossilPOST_PC_scores <- as.data.frame(TrunkFossilPOST_PCA)

POSTscores <-data.frame(TrunkPost.pca$x, species= (TrunkPost$species))

FossilPOST_PC_scores <-data.frame(TrunkFossilPOST_PCA, species= TrunkFossilT12$species)

All_POSTPC_scores <- (rbind(POSTscores, FossilPOST_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil

percentage_post <- round(TrunkPost.pca$sdev / sum(TrunkPost.pca$sdev) * 100, 2)# find percentage variance explained by PC's
percentage_post <- paste( colnames(POSTscores), "(", paste( as.character(percentage_post), "%", ")", sep="") )

# PLOT #
Post_plot<-ggplot(All_POSTPC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_post[1]) + ylab(percentage_post[2]) +
  scale_color_manual(name = "Species", breaks=levels(TrunkPost$species), values=c(speciescolors, fossilcolors)) + theme_classic() + ggtitle("Posterior Vertebrae")
Post_plot




# Sacral fossils

FossilSC_PCA <- predict(Sc.pca, TrunkFossilSc[,1:7])
FossilSC_PC_scores <- as.data.frame(FossilSC_PCA)

SCscores <- data.frame(Sc.pca$x, species= Sc$species)


FossilSC_PC_scores <- data.frame(FossilSC_PCA, species= TrunkFossilSc$species)

All_SCPC_scores <- (rbind(SCscores, FossilSC_PC_scores)) # create a new dataframe with the original PC scores and the PC scores of your fossil

percentage_Sc <- round(Sc.pca$sdev / sum(Sc.pca$sdev) * 100, 2)# find percentage variance explained by PC's
percentage_Sc <- paste( colnames(SCscores), "(", paste( as.character(percentage_Sc), "%", ")", sep="") )

# PLOT #
Sc_plot<-ggplot(All_SCPC_scores,aes(x=PC1,y=PC2,color=species)) + 
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(color=species), size = 1) +
  geom_point(size =2)+ xlab(percentage_Sc[1]) + ylab(percentage_Sc[2]) +
  scale_color_manual(name = "Species", breaks=levels(Sc$species), values=c(speciescolors, fossilcolors)) + theme_classic() + ggtitle("Sacral Vertebrae")
Sc_plot



# PLOT ALL TOGETHER #
par(oma=c(0,0,2,0))
gridExtra::grid.arrange(Ant_plot, Mid_plot ,Post_plot,Sc_plot,nrow = 2)





# T8 POSTERIOR EXTENSION MEASUREMENTS PLOTS #

library(EnvStats)
T8_extension

T8_extension<-mutate(T8_extension ,ratio = Cen_to_NeuAr/Cen_to_PoZy, .before = Cen_to_PoZy) #add new column for ratio

library(tidyr)
data_long <- gather(T8_extension, Type , Measurement, Cen_to_NeuAr:Cen_to_PoZy, factor_key=TRUE) #convert to long format

library(ggplot2)
s <- ggplot(data_long, aes(species, Measurement, fill = Type)) + geom_boxplot(position = "dodge") + theme_classic() + ylab("Value") +
  theme(legend.position="bottom") + scale_fill_discrete(name = "Measurement Type", labels= c("NSPE", "NSPE/PoZyPE", "PoZyPE"))
s





# Assess sample size per species
library(tidyverse)

PoAtVerts_wofossil %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(N = n())

# *removed A. subsalsum and A. ordinarium* 

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





### K NEAREST NEIGHBOR   ###:Non-parametric
library(caret)

#make KNN model using LOOCV to find optimal k

# SACRAL VERTEBRAE #
set.seed(123)
KNNmodel_sc <- train(
  species ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data = TrunkSc_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_sc) # plot accuracy vs k
KNNmodel_sc$bestTune # optimal k

predicted.classes_sc <- KNNmodel_sc %>% predict(TrunkSc_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_sc)
mean(predicted.classes_sc == TrunkSc_sub$species) #overall accuracy

# assess accuracy per species
accKNN_sc <- table(TrunkSc_sub$species,predicted.classes_sc)
accKNN_sc
diag(prop.table(accKNN_sc, 1))


# POSTERIOR VERTEBRAE #
set.seed(123)
KNNmodel_post <- train(
  species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data = TrunkPost_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_post) # plot accuracy vs k
KNNmodel_post$bestTune # optimal k

predicted.classes_post <- KNNmodel_post %>% predict(TrunkPost_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_post)
mean(predicted.classes_post == TrunkPost_sub$species) #overall accuracy

# assess accuracy per species
accKNN_post <- table(TrunkPost_sub$species,predicted.classes_post)
accKNN_post
diag(prop.table(accKNN_post, 1))


# MIDDLE VERTEBRAE #
set.seed(123)
KNNmodel_mid <- train(
  species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data = TrunkMid_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_mid) # plot accuracy vs k
KNNmodel_mid$bestTune # optimal k

predicted.classes_mid <- KNNmodel_mid %>% predict(TrunkMid_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_mid)
mean(predicted.classes_mid == TrunkMid_sub$species) #overall accuracy

# assess accuracy per species
accKNN_mid <- table(TrunkMid_sub$species,predicted.classes_mid)
accKNN_mid
diag(prop.table(accKNN_mid, 1))


# ANTERIOR VERTEBRAE #
set.seed(123)
KNNmodel_ant <- train(
  species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data = TrunkAnt_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_ant) # plot accuracy vs k
KNNmodel_ant$bestTune # optimal k

predicted.classes_ant <- KNNmodel_ant %>% predict(TrunkAnt_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_ant)
mean(predicted.classes_ant == TrunkAnt_sub$species) #overall accuracy

# assess accuracy per species
accKNN_ant <- table(TrunkAnt_sub$species,predicted.classes_ant)
accKNN_ant
diag(prop.table(accKNN_ant, 1))





## KNN FOSSIL CLASSIFICATIONS ##
library(class)

# ANTERIOR FOSSILS #

KnnAntPrediction <- knn(TrunkAnt_sub[,1:7], TrunkFossilAnt[,1:7],
                           TrunkAnt_sub$species, k=KNNmodel_ant$bestTune , prob=TRUE)
KnnAntPrediction

# t <- cbind(as.character(TrunkFossilAnt$species), as.character(KnnAntPrediction_k7), as.character(attr(KnnAntPrediction_k7, 'prob')))
# t <- as.data.frame(t)
# t$V3<- as.numeric(t$V3)
# t$V3 <- round(t$V3, 2)
# write.table(t, file = "Anterior fossils.txt", sep = ",", quote = FALSE, row.names = T)


# MIDDLE FOSSILS #
KnnMidPrediction <- knn(TrunkMid_sub[,1:7], TrunkFossilT8[,1:7],
                           TrunkMid_sub$species, k=KNNmodel_mid$bestTune, prob=TRUE)
KnnMidPrediction


# POSTERIOR FOSSILS #
KnnPostPrediction <- knn(TrunkPost_sub[,1:7], TrunkFossilT12[,1:7],
                            TrunkPost_sub$species, k=KNNmodel_post$bestTune, prob=TRUE)
KnnPostPrediction


# SACRAL VERTEBRAE FOSSILS #
KnnPostPrediction <- knn(TrunkSc_sub[,1:7], TrunkFossilSc[,1:7],
                            TrunkSc_sub$species, k=KNNmodel_sc$bestTune, prob=TRUE)
KnnPostPrediction






### RANDOM FOREST CLASSIFICATION ###:Non-parametric
library(randomForest)

# ANTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_ant <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkAnt_sub, importance=FALSE)
print(Atlas.rf_ant)
rf_acc_ant <- Atlas.rf_ant$confusion
rf_acc_ant <- 1-rf_acc_ant[,14] # percent correct classification
rf_acc_ant

mean(Atlas.rf_ant$predicted == TrunkAnt_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_ant = predict(Atlas.rf_ant, newdata = TrunkFossilAnt[,1:7])
y_pred_ant



# MIDDLE VERTEBRAE #
set.seed(123)
Atlas.rf_mid <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkMid_sub, importance=FALSE)
print(Atlas.rf_mid)
rf_acc_mid <- Atlas.rf_mid$confusion
rf_acc_mid <- 1-rf_acc_mid[,14] # percent correct classification
rf_acc_mid

mean(Atlas.rf_mid$predicted == TrunkMid_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_mid = predict(Atlas.rf_mid, newdata = TrunkFossilT8[,1:7])
y_pred_mid



# POSTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_post <- randomForest(species ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkPost_sub, importance=FALSE)
print(Atlas.rf_post)
rf_acc_post <- Atlas.rf_post$confusion
rf_acc_post <- 1-rf_acc_post[,14] # percent correct classification
rf_acc_post

mean(Atlas.rf_post$predicted == TrunkPost_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_post = predict(Atlas.rf_post, newdata = TrunkFossilT12[,1:7])
y_pred_post



# SACRAL VERTEBRAE #
set.seed(123)
Atlas.rf_sc <- randomForest(species ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data=TrunkSc_sub, importance=FALSE)
print(Atlas.rf_sc)
rf_acc_sc <- Atlas.rf_sc$confusion
rf_acc_sc <- 1-rf_acc_sc[,14] # percent correct classification
rf_acc_sc

mean(Atlas.rf_sc$predicted == TrunkSc_sub$species) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_sc = predict(Atlas.rf_sc, newdata = TrunkFossilSc[,1:7])
y_pred_sc






### MAHALAHOBIS DISTANCES AND NEIGHBOR JOINING ###
library(HDMD)
library("ape")
library(phytools)

Mah_Dis <- function(pairwiseMah) {
  names = rownames(pairwiseMah$means) #capture labels
  mahala = sqrt(pairwiseMah$distance) #mahalanobis distance
  rownames(mahala) = names #set rownames in the dissimilarity matrix
  colnames(mahala) = names #set colnames in the dissimilarity matrix
  return(mahala <- as.dist(mahala)) #this is the mahalanobis dissimilarity matrix 
} # return mahalanobis dissimilarity matrix


# ANTERIOR VERTEBRAE #
ANTTotal <- rbind(TrunkAnt[,1:8], TrunkFossilAnt)
Mahala1ANT = pairwise.mahalanobis(ANTTotal[,1:7], ANTTotal$species, digits = 3)
Mah_DisAnt <- Mah_Dis(Mahala1ANT)

trANT <- nj(Mah_DisAnt) #neighbor joining

plot((as.phylo(trANT)),type="unrooted",cex=0.6,
     use.edge.length=FALSE,lab4ut="axial",
     no.margin=TRUE)

# MIDDLE VERTEBRAE #
MIDTotal <- rbind(TrunkMid[,1:8], TrunkFossilT8)
Mahala1MID = pairwise.mahalanobis(MIDTotal[,1:7], MIDTotal$species, digits = 3)
Mah_DisMID <- Mah_Dis(Mahala1MID)

trMID <- nj(Mah_DisMID) #neighbor joining

plot((as.phylo(trMID)),type="unrooted",cex=0.6,
     use.edge.length=TRUE,lab4ut="axial",
     no.margin=TRUE)


# POSTERIOR VERTEBRAE #
POSTTotal <- rbind(TrunkPost[,1:8], TrunkFossilT12)
Mahala1POST = pairwise.mahalanobis(POSTTotal[,1:7], POSTTotal$species, digits = 3)
Mah_DisPOST <- Mah_Dis(Mahala1POST)

trPOST <- nj(Mah_DisPOST) #neighbor joining

plot((as.phylo(trPOST)),type="unrooted",cex=0.6,
     use.edge.length=TRUE,lab4ut="axial",
     no.margin=TRUE)# MIDDLE VERTEBRAE #


# SACRAL VERTEBRAE #
SCTotal <- rbind(Sc[,1:8], TrunkFossilSc)
Mahala1SC = pairwise.mahalanobis(SCTotal[,1:7], SCTotal$species, digits = 3)
Mah_DisSC <- Mah_Dis(Mahala1SC)

trSC <- nj(Mah_DisSC) #neighbor joining

plot((as.phylo(trSC)),type="unrooted",cex=0.6,
     use.edge.length=TRUE,lab4ut="axial",
     no.margin=TRUE)# MIDDLE VERTEBRAE #





# AMBYSTOMA CLADE CLASSIFICATION #

TrunkAnt_sub$clades <- recode(TrunkAnt_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
TrunkMid_sub$clades <- recode(TrunkMid_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
TrunkPost_sub$clades <- recode(TrunkPost_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")
TrunkSc_sub$clades <- recode(TrunkSc_sub$species, A.gracile = "A", A.talpoideum = "A", A.maculatum = "B", A.macrodactylum = "C", A.opacum = "D", A.laterale = "E", A.jeffersonianum = "E", A.mabeei = "F", A.texanum = "F", A.annulatum = "G", A.mavortium = "H", A.tigrinum = "H", A.velasci = "H")


### K NEAREST NEIGHBOR CLADES ###:Non-parametric
library(caret)

#make KNN model using LOOCV to find optimal k

# SACRAL VERTEBRAE #
set.seed(123)
KNNmodel_sc_clade <- train(
  clades ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data = TrunkSc_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_sc_clade) # plot accuracy vs k
KNNmodel_sc_clade$bestTune # optimal k

predicted.classes_sc_clade <- KNNmodel_sc_clade %>% predict(TrunkSc_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_sc_clade)
mean(predicted.classes_sc_clade == TrunkSc_sub$clades) #overall accuracy

# assess accuracy per clades
accKNN_sc_clade <- table(TrunkSc_sub$clades,predicted.classes_sc_clade)
accKNN_sc_clade
diag(prop.table(accKNN_sc_clade, 1))


# POSTERIOR VERTEBRAE #
set.seed(123)
KNNmodel_post_clade <- train(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data = TrunkPost_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_post_clade) # plot accuracy vs k
KNNmodel_post_clade$bestTune # optimal k

predicted.classes_post_clade <- KNNmodel_post_clade %>% predict(TrunkPost_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_post_clade)
mean(predicted.classes_post_clade == TrunkPost_sub$clades) #overall accuracy

# assess accuracy per clades
accKNN_post_clade <- table(TrunkPost_sub$clades,predicted.classes_post_clade)
accKNN_post_clade
diag(prop.table(accKNN_post_clade, 1))


# MIDDLE VERTEBRAE #
set.seed(123)
KNNmodel_mid_clade <- train(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data = TrunkMid_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_mid_clade) # plot accuracy vs k
KNNmodel_mid_clade$bestTune # optimal k

predicted.classes_mid_clade <- KNNmodel_mid_clade %>% predict(TrunkMid_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_mid_clade)
mean(predicted.classes_mid_clade == TrunkMid_sub$clades) #overall accuracy

# assess accuracy per clades
accKNN_mid_clade <- table(TrunkMid_sub$clades,predicted.classes_mid_clade)
accKNN_mid_clade
diag(prop.table(accKNN_mid_clade, 1))


# ANTERIOR VERTEBRAE #
set.seed(123)
KNNmodel_ant_clade <- train(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data = TrunkAnt_sub, method = "knn",
  trControl = trainControl("LOOCV", number =1),
  preProcess = c("center"), #center the data
  tuneLength = 15)

plot(KNNmodel_ant_clade) # plot accuracy vs k
KNNmodel_ant_clade$bestTune # optimal k

predicted.classes_ant_clade <- KNNmodel_ant_clade %>% predict(TrunkAnt_sub[,1:7]) # predict class based on KNN model
head(predicted.classes_ant_clade)
mean(predicted.classes_ant_clade == TrunkAnt_sub$clades) #overall accuracy

# assess accuracy per clades
accKNN_ant_clade <- table(TrunkAnt_sub$clades,predicted.classes_ant_clade)
accKNN_ant_clade
diag(prop.table(accKNN_ant_clade, 1))





## KNN FOSSIL CLADE CLASSIFICATIONS ##
library(class)

# ANTERIOR FOSSILS #

KnnAntPrediction <- knn(TrunkAnt_sub[,1:7], TrunkFossilAnt[,1:7],
                        TrunkAnt_sub$clades, k=KNNmodel_ant_clade$bestTune , prob=TRUE)
KnnAntPrediction

# t <- cbind(as.character(TrunkFossilAnt$clades), as.character(KnnAntPrediction_k7), as.character(attr(KnnAntPrediction_k7, 'prob')))
# t <- as.data.frame(t)
# t$V3<- as.numeric(t$V3)
# t$V3 <- round(t$V3, 2)
# write.table(t, file = "Anterior fossils.txt", sep = ",", quote = FALSE, row.names = T)


# MIDDLE FOSSILS #
KnnMidPrediction <- knn(TrunkMid_sub[,1:7], TrunkFossilT8[,1:7],
                        TrunkMid_sub$clades, k=KNNmodel_mid_clade$bestTune, prob=TRUE)
KnnMidPrediction


# POSTERIOR FOSSILS #
KnnPostPrediction <- knn(TrunkPost_sub[,1:7], TrunkFossilT12[,1:7],
                         TrunkPost_sub$clades, k=KNNmodel_post_clade$bestTune, prob=TRUE)
KnnPostPrediction


# SACRAL VERTEBRAE FOSSILS #
KnnPostPrediction <- knn(TrunkSc_sub[,1:7], TrunkFossilSc[,1:7],
                         TrunkSc_sub$clades, k=KNNmodel_sc_clade$bestTune, prob=TRUE)
KnnPostPrediction






### RANDOM FOREST CLADE CLASSIFICATION ###:Non-parametric
library(randomForest)

# ANTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_ant_clade <- randomForest(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkAnt_sub, importance=FALSE)
print(Atlas.rf_ant_clade)
rf_acc_ant_clade <- Atlas.rf_ant_clade$confusion
rf_acc_ant_clade <- 1-rf_acc_ant_clade[,14] # percent correct classification
rf_acc_ant_clade

mean(Atlas.rf_ant_clade$predicted == TrunkAnt_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_ant_clade = predict(Atlas.rf_ant_clade, newdata = TrunkFossilAnt[,1:7])
y_pred_ant_clade



# MIDDLE VERTEBRAE #
set.seed(123)
Atlas.rf_mid_clade <- randomForest(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkMid_sub, importance=FALSE)
print(Atlas.rf_mid_clade)
rf_acc_mid_clade <- Atlas.rf_mid_clade$confusion
rf_acc_mid_clade <- 1-rf_acc_mid_clade[,14] # percent correct classification
rf_acc_mid_clade

mean(Atlas.rf_mid_clade$predicted == TrunkMid_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_mid_clade = predict(Atlas.rf_mid_clade, newdata = TrunkFossilT8[,1:7])
y_pred_mid_clade



# POSTERIOR VERTEBRAE #
set.seed(123)
Atlas.rf_post_clade <- randomForest(clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkPost_sub, importance=FALSE)
print(Atlas.rf_post_clade)
rf_acc_post_clade <- Atlas.rf_post_clade$confusion
rf_acc_post_clade <- 1-rf_acc_post_clade[,14] # percent correct classification
rf_acc_post_clade

mean(Atlas.rf_post_clade$predicted == TrunkPost_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_post_clade = predict(Atlas.rf_post_clade, newdata = TrunkFossilT12[,1:7])
y_pred_post_clade



# SACRAL VERTEBRAE #
set.seed(123)
Atlas.rf_sc_clade <- randomForest(clades ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data=TrunkSc_sub, importance=FALSE)
print(Atlas.rf_sc_clade)
rf_acc_sc_clade <- Atlas.rf_sc_clade$confusion
rf_acc_sc_clade <- 1-rf_acc_sc_clade[,14] # percent correct classification
rf_acc_sc_clade

mean(Atlas.rf_sc_clade$predicted == TrunkSc_sub$clades) #overall accuracy

# FOSSIL CLASSIFICATION #
y_pred_sc_clade = predict(Atlas.rf_sc_clade, newdata = TrunkFossilSc[,1:7])
y_pred_sc_clade





# MEASUREMENT RELATIVE IMPORTANCE BASED ON RF #

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
} # function to plot measurement conditional permutation importance

# CONDITIONAL PERUMATATION IMPORTANCE # *long time to run

Atlas.rf_ant_imp <- cforest(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkAnt_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
create_crfplot(Atlas.rf_ant_imp, conditional = TRUE)

Atlas.rf_mid_imp <- cforest(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkMid_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
create_crfplot(Atlas.rf_mid_imp, conditional = TRUE)

Atlas.rf_post_imp <- cforest(
  clades ~ X1a + X2a + X3a + X4a + X5a + X6a + X7a, data=TrunkPost_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
create_crfplot(Atlas.rf_post_imp, conditional = TRUE)

Atlas.rf_sc_imp <- cforest(
  clades ~ X14 + X15 + X16 + X17 + X18 + X19 + X20, data=TrunkSc_sub,
  control = cforest_unbiased(mtry = 2, ntree = 500)
)
create_crfplot(Atlas.rf_sc_imp, conditional = TRUE)






