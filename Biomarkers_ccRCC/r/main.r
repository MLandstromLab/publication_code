#####################################
##  Landstrom Biomarker Project 1  ##
#####################################

## Packages
library(ggplot2)
library(PCAtools)
library(pheatmap)
library(ggforce) ##ggplot extension
library(randomForest)
library(dplyr)
library(gridExtra)
library(pROC)
library(corrplot)
library(caret)
library(WriteXLS)
library(STRINGdb)

################################################################################################################
################################################################################################################

#################################
##  Adjusted correlation plots ##
#################################

## Load data
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj.csv"                                 ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)            ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID"))]                                ## delete name column

preAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_preAdj.csv"                                 ## file path
preAdj_dat <- read.table(preAdj_fi, sep=",")                                                                    ## load to table (weird issue when using header=TRUE)
rownames(preAdj_dat) <- preAdj_dat[,1]                                                                          ## set rownames
preAdj_dat <- preAdj_dat[ , !(names(preAdj_dat) %in% c("V1"))]                                                  ## delete name column (v1 cos no header read)
names(preAdj_dat) <- names(ageAdj_dat)                                                                          ## copy names from ageAdj (checked to match)

## Target proteins/genes/whatever they are:
## Also a function to do the three pre/post/versus plots of each gene
targetGenes <- c("CD27", "CXL17", "CYR61", "EPHA2", "IGF1R", "ITGB5", "RET", "TGFR.2", "TNFRSF19", "TNFSF13", "WFDC2", "WISP.1")
ageCorrectionPlotter <- function(geneName){
  ## before correction plot
  preAdjNPX <- preAdj_dat[[geneName]]; preAdjAGE <- preAdj_dat[["Age"]]; preAdj_frame <- data.frame(PreAdjustedNPX = preAdjNPX, Age=preAdjAGE)
  ggplot(preAdj_frame, aes(x=Age, y=PreAdjustedNPX)) + geom_point() + ggtitle(paste(geneName, " pre-correction", sep=""))
  outPath1 <- paste("/home/alastairm/Desktop/Landstrom/results/AgeCorrectionPlots/", geneName, "_preCorrection.pdf", sep="")
  ggsave(outPath1)
  
  ## post correction plot
  ageAdjNPX <- ageAdj_dat[[geneName]]; ageAdjAGE <- ageAdj_dat[["Age"]]; ageAdj_frame <- data.frame(AgeAdjustedNPX = ageAdjNPX, Age=ageAdjAGE)
  ggplot(ageAdj_frame, aes(x=Age, y=AgeAdjustedNPX)) + geom_point() + ggtitle(paste(geneName, " age-corrected", sep=""))
  outPath2 <- paste("/home/alastairm/Desktop/Landstrom/results/AgeCorrectionPlots/", geneName, "_ageCorrection.pdf", sep="")
  ggsave(outPath2)
  
  ## pre VS post plot
  lowlim = 0; highlim = 0
  if(geneName=="CD27"){
    lowlim = 6; highlim = 10
  } else if(geneName=="CXL17"){
    lowlim = 1; highlim = 8
  } else if(geneName=="CYR61"){
    lowlim = 3; highlim = 10
  } else if(geneName=="EPHA2"){
    lowlim = 1; highlim = 5
  } else if(geneName=="IGF1R"){
    lowlim = 2; highlim = 6
  } else if(geneName=="ITGB5"){
    lowlim = 7; highlim = 10
  } else if(geneName=="RET"){
    lowlim = 4; highlim = 8
  } else if(geneName=="TGFR.2"){
    lowlim = 3; highlim = 8
  } else if(geneName=="TNFRSF19"){
    lowlim = 5; highlim = 11
  } else if(geneName=="TNFSF13"){
    lowlim = 6; highlim = 12
  } else if(geneName=="WFDC2"){
    lowlim = 5; highlim = 10
  } else if(geneName=="WISP.1"){
    lowlim = 4; highlim = 10
  }

  versus_frame <- data.frame(PreAdjustedNPX = preAdjNPX, AgeAdjustedNPX = ageAdjNPX)
  plotVersus <- ggplot(versus_frame, aes(x=PreAdjustedNPX, y=AgeAdjustedNPX)) + ylim(c(lowlim, highlim)) + xlim(c(lowlim,highlim)) + geom_point() + ggtitle(paste(geneName, " pre VS post-correction"))
  plotVersus <- plotVersus + coord_fixed()
  outPath3 <- paste("/home/alastairm/Desktop/Landstrom/results/AgeCorrectionPlots/", geneName, "_versusPlot.pdf", sep="")
  ggsave(outPath3)
}

## run the function on all the targets
sapply(targetGenes, FUN=ageCorrectionPlotter)

################################################################################################################
################################################################################################################

###################################
##  Protein-Protein interactions ##
###################################

## Anything above 70% correlation is considered highly correlated
## We want to extract these, then the associated P-Value for that interaction, and place in a table

coef_fi <- "/home/alastairm/Desktop/Landstrom/data/coef.csv"
coef_dat <- read.table(coef_fi, sep=';')

coefP_fi <- "/home/alastairm/Desktop/Landstrom/data/coef_p.csv"
coefP_dat <- read.table(coefP_fi, sep=";")

## Done using python
## project repo/python/proteinInteraction.py

################################################################################################################
################################################################################################################

#################
##  PCA        ##
#################

## load data
## reloading objects again just to start "clean" with each sub-task
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj.csv"                                 ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID"))]                                ## delete name column

## split into read-data and "phenotype" data
ageAdj_values <- ageAdj_dat[1:(length(ageAdj_dat)-14)]                                                          ## remove last 14 (phenotype)
ageAdj_values <- data.frame(t(ageAdj_values))                                                                   ## flip for PCAtools
ageAdj_phenotype <- ageAdj_dat[tail(names(ageAdj_dat), 14)]                                                     ## keep only last 14 (phenotype)
rownames(ageAdj_phenotype) <- names(ageAdj_values)
ageAdj_phenotype[,'PlateID']<-factor(ageAdj_phenotype[,'PlateID'])                                              ## need integer class labels as factors

p <- pca(ageAdj_values, metadata = ageAdj_phenotype, removeVar = 0.1)                                           ## pca YoooOOo
pdf(file="/home/alastairm/Desktop/Landstrom/results/DiseaseStatePCA.pdf")
biplot(p, lab = NULL, colby = 'group', title="Disease state PCA", hline = 0, vline = 0, legendPosition = 'right')
dev.off()

pdf(file="/home/alastairm/Desktop/Landstrom/results/PlatePCA.pdf")                                          ## pca YoooOOo
biplot(p, lab = NULL, colby = 'PlateID',
       title="Plate effect PCA", hline = 0, vline = 0, legendPosition = 'right')
dev.off()

################################################################################################################
################################################################################################################

########################
##  Heatmap           ##
########################

## Reloading data for cleanliness
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj.csv"                                 ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID"))]                                ## delete name column

## Metadata and Data values -- split from main frame
ageAdj_values <- ageAdj_dat[1:(length(ageAdj_dat)-14)]                                                          ## remove last 14 (phenotype)
ageAdj_values <- data.frame(t(ageAdj_values))                                                                   ## flip for PCAtools
ageAdj_phenotype <- ageAdj_dat[tail(names(ageAdj_dat), 14)]                                                     ## keep only last 14 (phenotype)
rownames(ageAdj_phenotype) <- names(ageAdj_values)
ageAdj_phenotype[,'PlateID']<-factor(ageAdj_phenotype[,'PlateID'])                                              ## need integer class labels as factors

## subsect phenotype data
sub_ageAdj_phenotype <- ageAdj_phenotype[, c("Age", "Sex", "group")]

## Row means, sort, subselect, re-add rownames
rowMeans <- data.frame(rowMeans(ageAdj_values))
sortedMeans <- rowMeans[order(-rowMeans$rowMeans.ageAdj_values.), , drop = FALSE]
top20_sorted <- data.frame(sortedMeans[1:20,]); temp <- head(rownames(sortedMeans),n=20); rownames(top20_sorted) <- temp
top50_sorted <- data.frame(sortedMeans[1:50,]); temp <- head(rownames(sortedMeans),n=50); rownames(top50_sorted) <- temp
top70_sorted <- data.frame(sortedMeans[1:70,]); temp <- head(rownames(sortedMeans),n=70); rownames(top70_sorted) <- temp

## subsect the dataframe rows with the requisite proteins
top20_values <- ageAdj_values[rownames(ageAdj_values) %in% rownames(top20_sorted),]
top50_values <- ageAdj_values[rownames(ageAdj_values) %in% rownames(top50_sorted),]
top70_values <- ageAdj_values[rownames(ageAdj_values) %in% rownames(top70_sorted),]

## Lists taken from wilcox test in customer data
top20_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF")
top50_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG")
top70_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG","CD160","ERBB3","ITGB5","FR.gamma","Gal.1","CPE","TCL1A","S100A4","ERBB4","hK8","FCRLB","PPY","TGFR.2","IGF1R","TRAIL","MAD.homolog.5","GPC1","MIC.A.B","WIF.1","GPNMB")
top20_values <- ageAdj_values[rownames(ageAdj_values) %in% top20_proteins,]
top50_values <- ageAdj_values[rownames(ageAdj_values) %in% top50_proteins,]
top70_values <- ageAdj_values[rownames(ageAdj_values) %in% top70_proteins,]

## plot the graph(s)
plot <- pheatmap(top20_values,
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
                 cellheight=5, cellwidth=5, fontsize=5, main='Significantly altered proteins: Top 20',
                 filename="/home/alastairm/Desktop/Landstrom/results/Top20Heatmap.pdf", width=20, height=20)
plot <- pheatmap(top50_values,
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
                 cellheight=5, cellwidth=5, fontsize=5, main="Significantly altered proteins: Top 50",
                 filename="/home/alastairm/Desktop/Landstrom/results/Top50Heatmap.pdf", width=20, height=20)
plot <- pheatmap(top70_values,
                 show_rownames = T, show_colnames = T, cluster_rows = T, cluster_cols = T, clustering_distance_cols = "euclidean",
                 trace = "none", density.info = "none", annotation_col = sub_ageAdj_phenotype, annotation_names_col = T,
                 cellheight=5, cellwidth=5, fontsize=5, main="Significantly altered proteins: Top 70",
                 filename="/home/alastairm/Desktop/Landstrom/results/Top70Heatmap.pdf", width=20, height=20)

################################################################################################################
################################################################################################################

#############################
##  Comparative proteomics ##
#############################

## Taking the statistically significant altered proteins and doing NPX value box plots of disease vs control
## NPX = age adjusted data readings
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj.csv"                                 ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID"))]                                ## delete name column

## split into disease and control
disease_values <- ageAdj_dat[1:134,]
control_values <- ageAdj_dat[135:245,]

## From customer data wilcox test, top 70 significant proteins
top70_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG","CD160","ERBB3","ITGB5","FR.gamma","Gal.1","CPE","TCL1A","S100A4","ERBB4","hK8","FCRLB","PPY","TGFR.2","IGF1R","TRAIL","MAD.homolog.5","GPC1","MIC.A.B","WIF.1","GPNMB")

## Boxlots done in python
## project repo/python/boxplotGenerator.py
#write.csv(disease_values, '/home/alastairm/Desktop/Landstrom/results/diseaseValues.csv')
#write.csv(control_values, '/home/alastairm/Desktop/Landstrom/results/controlValues.csv')

boxplot_fi <- "/home/alastairm/Desktop/Landstrom/results/boxplotData.csv"
boxplot_dat <- read.table(boxplot_fi, sep=",", header=FALSE)
names(boxplot_dat) <- c("Status","Protein","NPXValue")

## Plot facet of boxplots
for(i in 1:4){
  p <- ggplot(data = boxplot_dat, aes(x=Status, y=NPXValue)) 
  p <- p + geom_boxplot(aes(fill=Status))
  p <- p + geom_point(aes(y=NPXValue, group=Status), position = position_dodge(width=0.75))
  p <- p + facet_wrap_paginate( ~ Protein, scales="free", nrow=5, ncol=4, page=i)
  p <- p + xlab("Status") + ylab("NPX Value") + ggtitle("Significantly altered proteins: NPX value spread")
  ggsave(paste0("/home/alastairm/Desktop/Landstrom/results/boxplotPage", i, ".pdf"), width=15, height=15)
}

################################################################################################################
################################################################################################################

####################
##  Random Forest ##
####################

## Load data again
# ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj.csv"                                 ## file path
#ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
#rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames

ageAdj_dat <- readRDS(file.path("/home/tre/projects/landstrom/ageAdj_dat.rds"))

## RF limit so do on top 50, select columns + group column
top50_proteins <- c('group', "ESM1","MK","SYND1","FGFBP1","IL6","TFPI2","HGF","ANXA1","WFDC2","MetAP2","TNFRSF6B","RET","DLL1","CAIX",
                    "CXCL13","TNFRSF4","WISP1","FRalpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM8",
                    "VIM","AREG","GZMB","KLK13","hK11","VEGFR3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFNgammaR1","MSLN","CTSV",
                    "FURIN","MUC16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG")
ageAdj_top50 <- ageAdj_dat[ ,colnames(ageAdj_dat) %in% top50_proteins]
ageAdj_top50$group <- factor(ageAdj_top50$group)

## seed for reproducability
## Make model, basic plot (saved for modification)
set.seed(1337)
modelRF <- randomForest(group~., ageAdj_top50, ntree=100, importance=TRUE, nodesize=5, na.action=na.roughfix)
impData <- data.frame(modelRF$importance)
impData <- impData[order(impData$MeanDecreaseAccuracy, decreasing=TRUE),]
impData$Protein <- factor(rownames(impData))
rownames(impData) <- NULL

## Respecify Protein name as relation to value (so order is retained for MeanDecreaseAccuracy)
impData$Protein = factor(impData$Protein, levels=impData[order(impData$MeanDecreaseAccuracy), "Protein"])
MDAPlot <- ggplot(impData, aes(Protein, MeanDecreaseAccuracy)) +
  geom_point(stat='identity', size=3) + 
  coord_flip() +
  labs(title="Multivariate Analysis", subtitle="Random Forest, Disease VS Control") + theme_minimal()

## Respecify Protein name as relation to value (so order is retained for MeanDecreaseGini)
impData$Protein = factor(impData$Protein, levels=impData[order(impData$MeanDecreaseGini), "Protein"])
MDGPlot <- ggplot(impData, aes(Protein, MeanDecreaseGini)) + 
  geom_point(stat="identity", size=3) + 
  coord_flip() +
  labs(title="Multivariate Analysis", subtitle="Random Forest, Disease VS Control") + theme_minimal()

grid.arrange(MDAPlot, MDGPlot, nrow=1)
ggsave("/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/RandomForestPlot.pdf", arrangeGrob(MDAPlot, MDGPlot), width = 10, height = 9)


## Assess the performance of the model using cross-validation
# Split data to training and validation sets
set.seed(42)
model.data <- ageAdj_top50
model.data$class <- ifelse(model.data$group == 1, "D", "C")
model.data <- model.data[,colnames(model.data) != "group"] # Remove group column
splitSample <- createDataPartition(model.data$class, p = 0.8, list = FALSE)
training.data <- model.data[splitSample,]
validation.data <- model.data[-splitSample,]


# Run cross-validated RF
fitControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = T,
  classProbs = T)

set.seed(13)
rfFit <- train(class ~ ., data = training.data, 
                 method = "rf", 
                 trControl = fitControl)

# Predictions using the best model
rfPred <- predict(rfFit, newdata = validation.data, type = "prob")
roc <- roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
pdf("/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/RandomForestROC.pdf")
roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
dev.off()

roc.coords <- coords(roc, x = "best", ret = "all")
WriteXLS(roc.coords, "/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/RandomForest_ROC_metrics.xlsx")

################################################################################################################
################################################################################################################

####################################
##  Penalised logistic regression ##
####################################

## Load data again
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/modelData.csv"                                             ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep="\t", header=TRUE, check.names=FALSE)                                   ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat[,1] <- NULL

## Models
anxa1Model <- glm(group ~ ANXA1, data=ageAdj_dat, family="binomial")
caixModel <- glm(group ~ CAIX, data=ageAdj_dat, family="binomial")
egfModel <- glm(group ~ EGF, data=ageAdj_dat, family="binomial")
esm1Model <- glm(group ~ ESM1, data=ageAdj_dat, family="binomial")
fgfbp1Model <- glm(group ~ FGFBP1, data=ageAdj_dat, family="binomial")
il6Model <- glm(group ~ IL6, data=ageAdj_dat, family="binomial")
metap2Model <- glm(group ~ MetAP2, data=ageAdj_dat, family="binomial")
mkModel <- glm(group ~ MK, data=ageAdj_dat, family="binomial")
synd1Model <- glm(group ~ SYND1, data=ageAdj_dat, family="binomial")
tnfrsf6bModel <- glm(group ~ TNFRSF6B, data=ageAdj_dat, family="binomial")
txlnaModel <- glm(group ~ TXLNA, ageAdj_dat, family="binomial")
allModel <- glm(group ~ ., ageAdj_dat, family="binomial")

## Probabilities
anxa1Prob = predict(anxa1Model, newdata = ageAdj_dat, type = "response")
caixProb = predict(caixModel, newdata = ageAdj_dat, type = "response")
egfProb = predict(egfModel, newdata = ageAdj_dat, type = "response")
esm1Prob = predict(esm1Model, newdata = ageAdj_dat, type = "response")
fgfbp1Prob = predict(fgfbp1Model, newdata = ageAdj_dat, type = "response")
il6Prob = predict(anxa1Model, newdata = ageAdj_dat, type = "response")
metap2Prob = predict(metap2Model, newdata = ageAdj_dat, type = "response")
mkProb = predict(mkModel, newdata = ageAdj_dat, type = "response")
synd1Prob = predict(synd1Model, newdata = ageAdj_dat, type = "response")
tnfrsf6bProb = predict(tnfrsf6bModel, newdata = ageAdj_dat, type = "response")
txlnaProb = predict(txlnaModel, newdata = ageAdj_dat, type = "response")
allProb = predict(allModel, newdata=ageAdj_dat, type="response")

## ROCAUC Curves
anxa1ROC = multiclass.roc(ageAdj_dat$group ~ anxa1Prob, plot = TRUE, print.auc = FALSE, col="#e6194b")
caixROC = multiclass.roc(ageAdj_dat$group ~ caixProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#bcf60c')
egfROC = multiclass.roc(ageAdj_dat$group ~ egfProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#808000')
esm1ROC = multiclass.roc(ageAdj_dat$group ~ esm1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b')
fgfbp1ROC = multiclass.roc(ageAdj_dat$group ~ fgfbp1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe')
il6ROC = multiclass.roc(ageAdj_dat$group ~ il6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1')
metap2ROC = multiclass.roc(ageAdj_dat$group ~ metap2Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#469990')
mkROC = multiclass.roc(ageAdj_dat$group ~ mkProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#911eb4')
synd1ROC = multiclass.roc(ageAdj_dat$group ~ synd1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#800000')
tnfrsf6bROC = multiclass.roc(ageAdj_dat$group ~ tnfrsf6bProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#aaffc3')
txlnaROC = multiclass.roc(ageAdj_dat$group ~ txlnaProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#f58231')
allROC = multiclass.roc(ageAdj_dat$group ~ allProb, plot=TRUE, print.auc=FALSE, add=TRUE, col='#911eb4')

legend("bottomright", 
       legend = c("ANXA1", "CAIX", "EGF", "ESM1", "FGF.BP1", "IL6", "MetAP2", "MK", "SYND1", "TNFRSF6B", "TXLNA"), 
       col = c("#e6194b", "#bcf60c", "#808000", "#3cb44b", "#fabebe", '#ffd8b1', '#469990', '#911eb4', '#800000', '#aaffc3', '#f58231', '#911eb4'),
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

################################################################################################################
################################################################################################################

# The figure with the ROC curve currently has 11 different biomarkers that were increased in both tumor and control.  
# I have the following request:  
# Can you also make a ROC curve with  
# 1.) All the 11 (ANXA1, CAIX, EGF, ESM1, FGF.BP1, IL6, MetAP2, MK, SYND1, TNFRSF6B and TXLNA) biomarkers clubbed together,  
# 2.) 7 biomarkers (CAIX, ESM1, FGF.BP1, IL6, MK, SYND1 and TNFRSF6B) that were increased in tumors individually.  
# 3.) 7 biomarkers (CAIX, ESM1, FGF.BP1, IL6, MK, SYND1 and TNFRSF6B) that were increased in the tumors combined. 

# First read in probable plot data and try replicating the old ROC plot
ageAdj_dat <- read.csv("/home/data/project_data/landstrom_biomarker_project1/oncollAnalysis_ageAdj.csv", check.names = F)
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          
ageAdj_dat[,1] <- NULL
colnames(ageAdj_dat) <- gsub("\\.", "", colnames(ageAdj_dat))

# Convert group to binary vector 0/1
ageAdj_dat$group <- ifelse(ageAdj_dat$group == "Ctrl", 0, 1)
saveRDS(ageAdj_dat, file.path("/home/tre/projects/landstrom/ageAdj_dat.rds"))

# Subset to contain only columns of interest
ageAdj_dat <- ageAdj_dat[,c("ANXA1", "CAIX", "EGF", "ESM1", "FGFBP1", "IL6", "MetAP2", "MK", "SYND1", "TNFRSF6B", "TXLNA", "group")]


## Models
anxa1Model <- glm(group ~ ANXA1, data=ageAdj_dat, family="binomial")
caixModel <- glm(group ~ CAIX, data=ageAdj_dat, family="binomial")
egfModel <- glm(group ~ EGF, data=ageAdj_dat, family="binomial")
esm1Model <- glm(group ~ ESM1, data=ageAdj_dat, family="binomial")
fgfbp1Model <- glm(group ~ FGFBP1, data=ageAdj_dat, family="binomial")
il6Model <- glm(group ~ IL6, data=ageAdj_dat, family="binomial")
metap2Model <- glm(group ~ MetAP2, data=ageAdj_dat, family="binomial")
mkModel <- glm(group ~ MK, data=ageAdj_dat, family="binomial")
synd1Model <- glm(group ~ SYND1, data=ageAdj_dat, family="binomial")
tnfrsf6bModel <- glm(group ~ TNFRSF6B, data=ageAdj_dat, family="binomial")
txlnaModel <- glm(group ~ TXLNA, ageAdj_dat, family="binomial")
# allModel <- glm(group ~ ., ageAdj_dat, family="binomial")

## Probabilities
anxa1Prob = predict(anxa1Model, newdata = ageAdj_dat, type = "response")
caixProb = predict(caixModel, newdata = ageAdj_dat, type = "response")
egfProb = predict(egfModel, newdata = ageAdj_dat, type = "response")
esm1Prob = predict(esm1Model, newdata = ageAdj_dat, type = "response")
fgfbp1Prob = predict(fgfbp1Model, newdata = ageAdj_dat, type = "response")
il6Prob = predict(anxa1Model, newdata = ageAdj_dat, type = "response")
metap2Prob = predict(metap2Model, newdata = ageAdj_dat, type = "response")
mkProb = predict(mkModel, newdata = ageAdj_dat, type = "response")
synd1Prob = predict(synd1Model, newdata = ageAdj_dat, type = "response")
tnfrsf6bProb = predict(tnfrsf6bModel, newdata = ageAdj_dat, type = "response")
txlnaProb = predict(txlnaModel, newdata = ageAdj_dat, type = "response")
# allProb = predict(allModel, newdata=ageAdj_dat, type="response")

## ROCAUC Curves
anxa1ROC = roc(ageAdj_dat$group ~ anxa1Prob, plot = TRUE, print.auc = FALSE, col="#e6194b")
caixROC = roc(ageAdj_dat$group ~ caixProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#bcf60c')
egfROC = roc(ageAdj_dat$group ~ egfProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#808000')
esm1ROC = roc(ageAdj_dat$group ~ esm1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b')
fgfbp1ROC = roc(ageAdj_dat$group ~ fgfbp1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe')
il6ROC = roc(ageAdj_dat$group ~ il6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1')
metap2ROC = roc(ageAdj_dat$group ~ metap2Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#469990')
mkROC = roc(ageAdj_dat$group ~ mkProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#911eb4')
synd1ROC = roc(ageAdj_dat$group ~ synd1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#800000')
tnfrsf6bROC = roc(ageAdj_dat$group ~ tnfrsf6bProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#aaffc3')
txlnaROC = roc(ageAdj_dat$group ~ txlnaProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#f58231')
# allROC = roc(ageAdj_dat$group ~ allProb, plot=TRUE, print.auc=FALSE, add=TRUE, col='#911eb4')

legend(x = 0, y = 0.4, 
       legend = c("ANXA1", "CAIX", "EGF", "ESM1", "FGF.BP1", "IL6", "MetAP2", "MK", "SYND1", "TNFRSF6B", "TXLNA"), 
       col = c("#e6194b", "#bcf60c", "#808000", "#3cb44b", "#fabebe", '#ffd8b1', '#469990', '#911eb4', '#800000', '#aaffc3', '#f58231', '#911eb4'),
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
       cex = 0.7)


# 7 biomarkers that were increased
ageAdj_dat_7 <- ageAdj_dat[, c("CAIX", "ESM1", "FGFBP1", "IL6", "MK", "SYND1", "TNFRSF6B", "group")]

# Individually
caixROC = roc(ageAdj_dat_7$group ~ caixProb, plot = TRUE, print.auc = FALSE, col='#bcf60c')
esm1ROC = roc(ageAdj_dat_7$group ~ esm1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b')
fgfbp1ROC = roc(ageAdj_dat_7$group ~ fgfbp1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe')
il6ROC = roc(ageAdj_dat_7$group ~ il6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1')
mkROC = roc(ageAdj_dat_7$group ~ mkProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#911eb4')
synd1ROC = roc(ageAdj_dat_7$group ~ synd1Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#800000')
tnfrsf6bROC = roc(ageAdj_dat_7$group ~ tnfrsf6bProb, plot = TRUE, print.auc = FALSE, add=TRUE, col='#aaffc3')

legend(x = 0, y = 0.25,  
       legend = c("CAIX", "ESM1", "FGF.BP1", "IL6", "MK", "SYND1", "TNFRSF6B"), 
       col = c("#bcf60c", "#3cb44b", "#fabebe", '#ffd8b1', '#911eb4', '#800000', '#aaffc3'),
       lty = c(1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1),
       cex = 0.7)

# All 7 combined in one ROC
allModel7 <- glm(group ~ ., ageAdj_dat_7, family="binomial")
allProb7 = predict(allModel7, newdata = ageAdj_dat_7, type="response")

all7ROC = roc(ageAdj_dat_7$group ~ allProb7, plot = TRUE, print.auc = FALSE, col='red')


####################################
##  Penalised Logistic regression ##
####################################

library("pROC")
library("glmnet")
library("caret")
library(WriteXLS)

# Read in the data object stored above
ageAdj_dat <- readRDS(file.path("/home/tre/projects/landstrom/ageAdj_dat.rds"))

# Subset to contain only columns of interest, i.e. the previously identified best markers (performed by Felipe)
#ageAdj_dat <- ageAdj_dat[,c("ANXA1", "CAIX", "EGF", "ESM1", "FGFBP1", "IL6", "MetAP2", "MK", "SYND1", "TNFRSF6B", "TXLNA", "group")]
ageAdj_dat <- ageAdj_dat[,c("ANXA1", "EGF", "ESM1", "FGFBP1", "MetAP2", "MK", "SYND1", "group")]
all_glm <- glmnet(as.matrix(ageAdj_dat), ageAdj_dat$group, family="binomial")
plot(all_glm, xvar='lambda', label=TRUE)

# Split data for glmmet
set.seed(135337)

splitSample <- createDataPartition(ageAdj_dat$group, p = 0.8, list = FALSE)
training_expression <- ageAdj_dat[splitSample,]
training_phenotype <- ageAdj_dat$group[splitSample]
validation_expression <- ageAdj_dat[-splitSample,]
validation_phenotype <- data.frame(Sample = rownames(validation_expression), group = validation_expression$group)

# Run cross-validated glmnet using previously defined penalisation proportion alpha
# elasticnet <- cv.glmnet(as.matrix(training_expression[,-12]), training_expression$group, 
#                        family="binomial", nfolds=10, alpha=0.2)
elasticnet <- cv.glmnet(as.matrix(training_expression[,-8]), training_expression$group, 
                        family="binomial", nfolds=10, alpha=0.2)
plot(elasticnet)
# predict_validation <- predict(elasticnet, newx = as.matrix(validation_expression[,-12]), 
#                              s = c(elasticnet$lambda.min), type = "class")
predict_validation <- predict(elasticnet, newx = as.matrix(validation_expression[,-8]), 
                              s = c(elasticnet$lambda.min), type = "class")

validation_group <- as.matrix(validation_phenotype[,2])
row.names(validation_group) <- row.names(predict_validation)
colnames(validation_group) <- colnames(predict_validation)

# Get regression coefficients and order by decreasing absolute coef
coefs <- coef(elasticnet)
coefs <- data.frame(coefs[-1,]) # Remove intercept
coefs <- coefs[order(abs(coefs[,1]), decreasing = T),, drop = F]

# Make ROCs by starting from the protein with largest coefficient and adding one protein at time
ord.markers <- rownames(coefs)

## Models using full data set - 11 markers
# anxa1Model <- glm(group ~ ANXA1, data = ageAdj_dat, family="binomial")
# comb2Model <- glm(group ~ ANXA1 + ESM1, data = ageAdj_dat, family="binomial")
# comb3Model <- glm(group ~ ANXA1 + ESM1 + MetAP2, data = ageAdj_dat, family="binomial")
# comb4Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1, data = ageAdj_dat, family="binomial")
# comb5Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK, data = ageAdj_dat, family="binomial")
# comb6Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK + SYND1, data = ageAdj_dat, family="binomial")
# comb7Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK + SYND1 + EGF, data = ageAdj_dat, family="binomial")
# comb8Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK + SYND1 + EGF + TNFRSF6B, data = ageAdj_dat, family="binomial")
# comb9Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK + SYND1 + EGF + TNFRSF6B + TXLNA, data = ageAdj_dat, family="binomial")
# comb10Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK + SYND1 + EGF + TNFRSF6B + TXLNA + IL6, data = ageAdj_dat, family="binomial")
# comb11Model <- glm(group ~ ANXA1 + ESM1 + MetAP2 + FGFBP1 + MK + SYND1 + EGF + TNFRSF6B + TXLNA + IL6 + CAIX, data = ageAdj_dat, family="binomial")

## Models using full data set - 7 markers
anxa1Model <- glm(group ~ ANXA1, data = ageAdj_dat, family="binomial")
comb2Model <- glm(group ~ ANXA1 + ESM1, data = ageAdj_dat, family="binomial")
comb3Model <- glm(group ~ ANXA1 + ESM1 + FGFBP1, data = ageAdj_dat, family="binomial")
comb4Model <- glm(group ~ ANXA1 + ESM1 + FGFBP1 + SYND1, data = ageAdj_dat, family="binomial")
comb5Model <- glm(group ~ ANXA1 + ESM1 + FGFBP1 + SYND1 + MetAP2, data = ageAdj_dat, family="binomial")
comb6Model <- glm(group ~ ANXA1 + ESM1 + FGFBP1 + SYND1 + MetAP2 + MK, data = ageAdj_dat, family="binomial")
comb7Model <- glm(group ~ ANXA1 + ESM1 + FGFBP1 + SYND1 + MetAP2 + MK + EGF, data = ageAdj_dat, family="binomial")


# Probabilities
anxa1Prob = predict(anxa1Model, newdata = ageAdj_dat, type = "response")
comb2Prob = predict(comb2Model, newdata = ageAdj_dat, type = "response")
comb3Prob = predict(comb3Model, newdata = ageAdj_dat, type = "response")
comb4Prob = predict(comb4Model, newdata = ageAdj_dat, type = "response")
comb5Prob = predict(comb5Model, newdata = ageAdj_dat, type = "response")
comb6Prob = predict(comb6Model, newdata = ageAdj_dat, type = "response")
comb7Prob = predict(comb7Model, newdata = ageAdj_dat, type = "response")
# comb8Prob = predict(comb8Model, newdata = ageAdj_dat, type = "response")
# comb9Prob = predict(comb9Model, newdata = ageAdj_dat, type = "response")
# comb10Prob = predict(comb10Model, newdata = ageAdj_dat, type = "response")
# comb11Prob = predict(comb11Model, newdata = ageAdj_dat, type = "response")

# ROCAUC Curves
anxa1ROC = roc(ageAdj_dat$group ~ anxa1Prob, plot = TRUE, print.auc = FALSE, col="#e6194b")
comb2ROC = roc(ageAdj_dat$group ~ comb2Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#bcf60c')
comb3ROC = roc(ageAdj_dat$group ~ comb3Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#808000')
comb4ROC = roc(ageAdj_dat$group ~ comb4Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#3cb44b')
comb5ROC = roc(ageAdj_dat$group ~ comb5Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#fabebe')
comb6ROC = roc(ageAdj_dat$group ~ comb6Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#ffd8b1')
comb7ROC = roc(ageAdj_dat$group ~ comb7Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#469990')
# comb8ROC = roc(ageAdj_dat$group ~ comb8Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#911eb4')
# comb9ROC = roc(ageAdj_dat$group ~ comb9Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#800000')
# comb10ROC = roc(ageAdj_dat$group ~ comb10Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#aaffc3')
# comb11ROC = roc(ageAdj_dat$group ~ comb11Prob, plot = TRUE, print.auc = FALSE, add=TRUE, col='#f58231')

legend(x = 0, y = 0.5,  
       # legend = c("ANXA1", "comb2", "comb3", "comb4", "comb5", "comb6", "comb7", "comb8", "comb9", "comb10", "comb11"), 
       legend = c("ANXA1", "comb2", "comb3", "comb4", "comb5", "comb6", "comb7"),
       col = c("#e6194b", "#bcf60c", "#808000", "#3cb44b", "#fabebe", '#ffd8b1', '#469990', '#911eb4', '#800000', '#aaffc3', '#f58231', '#911eb4'),
       lty = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
       lwd = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
       cex = 0.7)

# Prepare result table - 11 markers
# res.table <- data.frame("Protein" = rownames(coefs),
#                         "Coef" = coefs[,1],
#                         "Combined numbers" = c("/", paste0("comb", seq(2,11))),
#                         "AUC" = NA,
#                         "Sen" = NA,
#                         "Spe" = NA,
#                         "ROC test P-value" = NA)
# 
# prob.list <- list(anxa1Prob, comb2Prob, comb3Prob, comb4Prob, comb5Prob, comb6Prob, comb7Prob, comb8Prob, comb9Prob, comb10Prob, comb11Prob)
# roc.list <- list(anxa1ROC, comb2ROC, comb3ROC, comb4ROC, comb5ROC, comb6ROC, comb7ROC, comb8ROC, comb9ROC, comb10ROC, comb11ROC)

# Prepare result table - 7 markers
res.table <- data.frame("Protein" = rownames(coefs),
                        "Coef" = coefs[,1],
                        "Combined numbers" = c("/", paste0("comb", seq(2,7))),
                        "AUC" = NA,
                        "Sen" = NA,
                        "Spe" = NA,
                        "ROC test P-value" = NA)

prob.list <- list(anxa1Prob, comb2Prob, comb3Prob, comb4Prob, comb5Prob, comb6Prob, comb7Prob)
roc.list <- list(anxa1ROC, comb2ROC, comb3ROC, comb4ROC, comb5ROC, comb6ROC, comb7ROC)

# AUC and other metrics for the first marker gene
AUC <- as.numeric(anxa1ROC$auc)
anxa1Pred <- ifelse(anxa1Prob > 0.5, "1", "0")
anxa1CM <- confusionMatrix(table(anxa1Pred, ageAdj_dat$group))
Sen <- anxa1CM$byClass["Sensitivity"]
Spe <- anxa1CM$byClass["Specificity"]

res.table[1,"AUC"] <- AUC
res.table[1,"Sen"] <- Sen
res.table[1,"Spe"] <- Spe

# Add the remaining values using a loop
# for(i in 2:11) {
for(i in 2:7) {  
  res.table[i,"AUC"] <- as.numeric(roc.list[[i]]$auc)
  pred <- ifelse(prob.list[[i]] > 0.5, "1", "0")
  cm <- confusionMatrix(table(pred, ageAdj_dat$group))
  res.table[i,"Sen"] <- cm$byClass["Sensitivity"]
  res.table[i,"Spe"] <- cm$byClass["Specificity"]
  
}

# Add ROC test p-values
res.table[1, "ROC.test.P.value"] <- "/"
res.table[2, "ROC.test.P.value"] <- "/"
res.table[3, "ROC.test.P.value"] <- roc.test(roc.list[[2]], roc.list[[3]])$p.value
res.table[4, "ROC.test.P.value"] <- roc.test(roc.list[[3]], roc.list[[4]])$p.value
res.table[5, "ROC.test.P.value"] <- roc.test(roc.list[[4]], roc.list[[5]])$p.value
res.table[6, "ROC.test.P.value"] <- roc.test(roc.list[[5]], roc.list[[6]])$p.value
res.table[7, "ROC.test.P.value"] <- roc.test(roc.list[[6]], roc.list[[7]])$p.value
# res.table[8, "ROC.test.P.value"] <- roc.test(roc.list[[7]], roc.list[[8]])$p.value
# res.table[9, "ROC.test.P.value"] <- roc.test(roc.list[[8]], roc.list[[9]])$p.value
# res.table[10, "ROC.test.P.value"] <- roc.test(roc.list[[9]], roc.list[[10]])$p.value
# res.table[11, "ROC.test.P.value"] <- roc.test(roc.list[[10]], roc.list[[11]])$p.value

# WriteXLS(res.table, file.path("/home/tre/projects/landstrom/combined_markers_results.xlsx"))
WriteXLS(res.table, file.path("/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/combined_markers_results_for_7.xlsx"))

## Assess the performance of the model with 7 biomarkers using cross-validation as was done above for the top 50 significant proteins
# Split data to training and validation sets
set.seed(42)
model.data <- ageAdj_dat
model.data$class <- ifelse(model.data$group == 1, "D", "C")
model.data <- model.data[,colnames(model.data) != "group"] # Remove group column
splitSample <- createDataPartition(model.data$class, p = 0.8, list = FALSE)
training.data <- model.data[splitSample,]
validation.data <- model.data[-splitSample,]

# Run cross-validated RF
fitControl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = T,
  classProbs = T)

set.seed(13)
rfFit <- train(class ~ ., data = training.data, 
               method = "rf", 
               trControl = fitControl)

# Predictions using the best model
rfPred <- predict(rfFit, newdata = validation.data, type = "prob")
roc <- roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
pdf("/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/RandomForestROC_for_7_biomarkers.pdf")
roc(validation.data$class ~ rfPred$D, plot = TRUE, print.auc = TRUE)
dev.off()

roc.coords <- coords(roc, x = "best", ret = "all")
WriteXLS(roc.coords, "/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/RandomForest_ROC_metrics_for_7_biomarkers.xlsx")


#####################################################################
## Modifications to previously generated result tables -  by Reija ##
#####################################################################

library(readr)
library(readxl)
library(tidyverse)
library(WriteXLS)

# Read in top70 table in which the p-value column had been modified to have . instead , as decimal separator 
top70signif_withStats <- read_excel("/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/top70signif_withStats_NEW.xlsx")

# Add group means and SDs
idx.D <- grep("^D_", colnames(top70signif_withStats))
idx.C <- grep("^ctrl", colnames(top70signif_withStats))

means.D <- rowMeans(top70signif_withStats[,idx.D])
sd.D <- sapply(1:70, function(x) sd(top70signif_withStats[x,idx.D]))
means.C <- rowMeans(top70signif_withStats[,idx.C])
sd.C <- sapply(1:70, function(x) sd(top70signif_withStats[x,idx.C]))
top70signif_withStats <- add_column(top70signif_withStats, "Mean of disease samples" = means.D, .after = "Standard deviation")
top70signif_withStats <- add_column(top70signif_withStats, "SD of disease samples" = sd.D, .after = "Mean of disease samples")
top70signif_withStats <- add_column(top70signif_withStats, "Mean of control samples" = means.C, .after = "SD of disease samples")
top70signif_withStats <- add_column(top70signif_withStats, "SD of control samples" = sd.C, .after = "Mean of control samples")

WriteXLS(top70signif_withStats, "/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/top70signif_withGroupStats.xlsx", BoldHeaderRow = T)
write_csv(top70signif_withStats, "/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/top70signif_withGroupStats.csv")


####################################################
## Checking STRINGdb for interactions -  by Reija ##
####################################################

int.prot <- colnames(ageAdj_dat)
int.prot <- head(int.prot, -1) # Remove group column

# Add current gene symbols to some proteins
prot.df <- data.frame("Name" = int.prot,
                      "Symbol" = int.prot)
prot.df$Symbol[prot.df$Symbol == "MetAP2"]  <- "METAP2"
prot.df$Symbol[prot.df$Symbol == "MK"]  <- "MDK"

# Load database and map to String IDs
string_db <- STRINGdb$new(version = "11.5", species = 9606, network_type = "full", input_directory = "")
mapping <- string_db$map(prot.df, "Symbol")

interactions <- string_db$get_interactions(mapping$STRING_id)

interactions$from_Symbol <- mapping$Symbol[match(interactions$from, mapping$STRING_id)]
interactions$to_Symbol <- mapping$Symbol[match(interactions$to, mapping$STRING_id)]

string_db$plot_network(mapping$STRING_id)

# Read in full database
StringDB <- read.table("/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/9606.protein.links.full.v11.5.txt", header = T)

StringDB.exp <- filter(StringDB, experiments != 0)
StringDB.exp.filt <- filter(StringDB.exp, protein1 %in% mapping$STRING_id | protein2 %in% mapping$STRING_id)

StringDB.exp.filt$protein1_Symbol <- mapping$Symbol[match(StringDB.exp.filt$protein1, mapping$STRING_id)]
StringDB.exp.filt$protein2_Symbol <- mapping$Symbol[match(StringDB.exp.filt$protein2, mapping$STRING_id)]

WriteXLS(StringDB.exp.filt, "/home/data/project_results/landstrom_biomarker_project1/latest_requests_2023/StringDB_exp_interactions.xlsx",
         BoldHeaderRow = T)

############################
##  Clinical Correlation  ##
############################

## Load data again
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_clinical.csv"                               ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
rownames(ageAdj_dat) <- ageAdj_dat$CaseID                                                                       ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "CaseID"))]                                  ## delete name column

## Load clinical data
clinical_fi <- "/home/alastairm/Desktop/Landstrom/data/clinicalCSV.csv"
clinical_dat <- read.table(clinical_fi, sep=",", header=TRUE, check.names=FALSE)

## Top 20/50/70 protein markers
## whoops ignore the above -- lists taken from wilcox test in customer data
top20_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF")
top50_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG")
top70_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG","CD160","ERBB3","ITGB5","FR.gamma","Gal.1","CPE","TCL1A","S100A4","ERBB4","hK8","FCRLB","PPY","TGFR.2","IGF1R","TRAIL","MAD.homolog.5","GPC1","MIC.A.B","WIF.1","GPNMB")
top20_values <- subset(ageAdj_dat, select=top20_proteins)
top50_values <- subset(ageAdj_dat, select=top50_proteins)
top70_values <- subset(ageAdj_dat, select=top70_proteins)

## plot colours
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

## Match ageAdjusted with clinical BY caseID
top20_mergedClinical <- merge(top20_values, clinical_dat, by.x = "row.names", by.y = "CaseID")
rownames(top20_mergedClinical) <- top20_mergedClinical[,1]; top20_mergedClinical[,1] <- NULL
top50_mergedClinical <- merge(top50_values, clinical_dat, by.x = "row.names", by.y = "CaseID")
rownames(top50_mergedClinical) <- top50_mergedClinical[,1]; top50_mergedClinical[,1] <- NULL
top70_mergedClinical <- merge(top70_values, clinical_dat, by.x = "row.names", by.y = "CaseID")
rownames(top70_mergedClinical) <- top70_mergedClinical[,1]; top70_mergedClinical[,1] <- NULL

## correlation matrices
top20_corrMatrix <- cor(top20_mergedClinical)
top50_corrMatrix <- cor(top50_mergedClinical)
top70_corrMatrix <- cor(top70_mergedClinical)
## save to CSV
write.csv(top20_corrMatrix, '/home/alastairm/Desktop/Landstrom/results/correlationPlotsTweaked/top20_corrMatrix.csv')
write.csv(top50_corrMatrix, '/home/alastairm/Desktop/Landstrom/results/correlationPlotsTweaked/top50_corrMatrix.csv')
write.csv(top70_corrMatrix, '/home/alastairm/Desktop/Landstrom/results/correlationPlotsTweaked/top70_corrMatrix.csv')

## Raw plots
top20_corrPlot <- corrplot(top20_corrMatrix, type="lower")
top50_corrPlot <- corrplot(top50_corrMatrix, type="lower")
top70_corrPlot <- corrplot(top70_corrMatrix, type="lower")

## Plots with hclust ordering
top20_corrPlotHclust <- corrplot(top20_corrMatrix, type="lower", order="hclust")
top50_corrPlotHclust <- corrplot(top50_corrMatrix, type="lower", order="hclust")
top70_corrPlotHclust <- corrplot(top70_corrMatrix, type="lower", order="hclust")

## Function for mapping p-values to matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

## corr + signif matrices
top20_significance <- cor.mtest(top20_mergedClinical)
top50_significance <- cor.mtest(top50_mergedClinical)
top70_significance <- cor.mtest(top70_mergedClinical)
## save to CSV
write.csv(top20_significance, '/home/alastairm/Desktop/Landstrom/results/correlationPlotsTweaked/top20_corrSignificance.csv')
write.csv(top50_significance, '/home/alastairm/Desktop/Landstrom/results/correlationPlotsTweaked/top50_corrSignificance.csv')
write.csv(top70_significance, '/home/alastairm/Desktop/Landstrom/results/correlationPlotsTweaked/top70_corrSignificance.csv')

## Plot with all info
top20_signifHclust <- corrplot(top20_corrMatrix, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)
top50_signifHclust <- corrplot(top50_corrMatrix, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)
top70_signifHclust <- corrplot(top70_corrMatrix, col=col(200), order="hclust",
                               tl.col="black", tl.srt=45, diag=FALSE)

################################################################################################################
################################################################################################################

####################################
##  Table of significant proteins ##
####################################

## Reloading data for cleanliness
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj.csv"                                 ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "sampleID"))]                                ## delete name column

## Metadata and Data values -- split from main frame
ageAdj_values <- ageAdj_dat[1:(length(ageAdj_dat)-14)]                                                          ## remove last 14 (phenotype)
ageAdj_values <- data.frame(t(ageAdj_values))                                                                   ## flip for PCAtools
ageAdj_phenotype <- ageAdj_dat[tail(names(ageAdj_dat), 14)]                                                     ## keep only last 14 (phenotype)
rownames(ageAdj_phenotype) <- names(ageAdj_values)
ageAdj_phenotype[,'PlateID']<-factor(ageAdj_phenotype[,'PlateID'])                                              ## need integer class labels as factors

## statistical data
stats_fi <- "/home/alastairm/Desktop/Landstrom/Customer/Material\ from\ customer/wilcoxTest.csv"
stats_dat <- read.table(stats_fi, sep=";", header=TRUE, check.names=FALSE)
colnames(stats_dat)[1] <- "protein"

## Lists taken from wilcox test in customer data
top20_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF")
top50_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG")
top70_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG","CD160","ERBB3","ITGB5","FR.gamma","Gal.1","CPE","TCL1A","S100A4","ERBB4","hK8","FCRLB","PPY","TGFR.2","IGF1R","TRAIL","MAD.homolog.5","GPC1","MIC.A.B","WIF.1","GPNMB")
top20_values <- ageAdj_values[rownames(ageAdj_values) %in% top20_proteins,]
top50_values <- ageAdj_values[rownames(ageAdj_values) %in% top50_proteins,]
top70_values <- ageAdj_values[rownames(ageAdj_values) %in% top70_proteins,]

## calculate mean value for each protein (row)
top20_values$meanVal <- rowMeans(top20_values, na.rm=TRUE)
top50_values$meanVal <- rowMeans(top50_values, na.rm=TRUE)
top70_values$meanVal <- rowMeans(top70_values, na.rm=TRUE)

## calculate sd for each row
top20_values$std.dev <- apply(top20_values,1,sd)
top50_values$std.dev <- apply(top50_values,1,sd)
top70_values$std.dev <- apply(top70_values,1,sd)

## grab p.val for each protein and plop it in
stats_top20sub <- subset(stats_dat, stats_dat$protein %in% top20_proteins)
stats_top20sub <- stats_top20sub[order(match(stats_top20sub$protein, rownames(top20_values))),]
top20_values$p.val <- stats_top20sub$p

stats_top50sub <- subset(stats_dat, stats_dat$protein %in% top50_proteins)
stats_top50sub <- stats_top50sub[order(match(stats_top50sub$protein, rownames(top50_values))),]
top50_values$p.val <- stats_top50sub$p

stats_top70sub <- subset(stats_dat, stats_dat$protein %in% top70_proteins)
stats_top70sub <- stats_top70sub[order(match(stats_top70sub$protein, rownames(top70_values))),]
top70_values$p.val <- stats_top70sub$p

## write to file
write.csv(top20_values, '/home/alastairm/Desktop/Landstrom/results/top20signif_withStats.csv')
write.csv(top50_values, '/home/alastairm/Desktop/Landstrom/results/top50signif_withStats.csv')
write.csv(top70_values, '/home/alastairm/Desktop/Landstrom/results/top70signif_withStats.csv')

################################################################################################################
################################################################################################################

###########################################################
##  Table of associated proteins with clinical variables ##
###########################################################

## Load data again
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_clinical.csv"                               ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                    ## load to table
rownames(ageAdj_dat) <- ageAdj_dat$CaseID                                                                       ## set rownames
ageAdj_dat <- ageAdj_dat[ , !(names(ageAdj_dat) %in% c("SampleID", "CaseID"))]                                  ## delete name column

## Load clinical data
clinical_fi <- "/home/alastairm/Desktop/Landstrom/data/clinicalCSV.csv"
clinical_dat <- read.table(clinical_fi, sep=",", header=TRUE, check.names=FALSE)

## Top 20/50/70 protein markers
## Lists taken from wilcox test in customer data
top20_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF")
top50_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG")
top70_proteins <- c("ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG","CD160","ERBB3","ITGB5","FR.gamma","Gal.1","CPE","TCL1A","S100A4","ERBB4","hK8","FCRLB","PPY","TGFR.2","IGF1R","TRAIL","MAD.homolog.5","GPC1","MIC.A.B","WIF.1","GPNMB")
top20_values <- subset(ageAdj_dat, select=top20_proteins)
top50_values <- subset(ageAdj_dat, select=top50_proteins)
top70_values <- subset(ageAdj_dat, select=top70_proteins)

## Match ageAdjusted with clinical BY caseID
top20_mergedClinical <- merge(top20_values, clinical_dat, by.x = "row.names", by.y = "CaseID")
rownames(top20_mergedClinical) <- top20_mergedClinical[,1]; top20_mergedClinical[,1] <- NULL
top50_mergedClinical <- merge(top50_values, clinical_dat, by.x = "row.names", by.y = "CaseID")
rownames(top50_mergedClinical) <- top50_mergedClinical[,1]; top50_mergedClinical[,1] <- NULL
top70_mergedClinical <- merge(top70_values, clinical_dat, by.x = "row.names", by.y = "CaseID")
rownames(top70_mergedClinical) <- top70_mergedClinical[,1]; top70_mergedClinical[,1] <- NULL

## correlation matrices
#test <- cor.test(top20_mergedClinical$ESM.1, top20_mergedClinical$Sex, method="pearson")

## Function for mapping p-values to matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

## corr + signif matrices
top20_significance <- cor.mtest(top20_mergedClinical)
top50_significance <- cor.mtest(top50_mergedClinical)
top70_significance <- cor.mtest(top70_mergedClinical)

## Function to remove non-significant values
remove_95 = function(x){
  ifelse(x > 0.05, NA, x)
}

## removed non signif values
top20_prunedSignif <- apply(top20_significance, c(1,2), FUN = remove_95)
top50_prunedSignif <- apply(top50_significance, c(1,2), FUN = remove_95)
top70_prunedSignif <- apply(top70_significance, c(1,2), FUN = remove_95)

## Write.csv
write.csv(top20_prunedSignif, '/home/alastairm/Desktop/Landstrom/results/top20associations_signif.csv')
write.csv(top50_prunedSignif, '/home/alastairm/Desktop/Landstrom/results/top50associations_signif.csv')
write.csv(top70_prunedSignif, '/home/alastairm/Desktop/Landstrom/results/top70associations_signif.csv')

## get protein names per clinical value (row)
top20_signifAssocProteins <- apply(top20_prunedSignif,1,function(x) names(which(!is.na(x))))
top50_signifAssocProteins <- apply(top50_prunedSignif,1,function(x) names(which(!is.na(x))))
top70_signifAssocProteins <- apply(top70_prunedSignif,1,function(x) names(which(!is.na(x))))

################################################################################################################
################################################################################################################

########################################
##  Volcano plots for top 70 proteins ##
########################################

## Load data again
ageAdj_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj_DDS.csv"                             ## file path
ageAdj_dat <- read.table(ageAdj_fi, sep=",", header=TRUE, check.names=FALSE)                                   ## load to table
rownames(ageAdj_dat) <- ageAdj_dat[,1]                                                                          ## set rownames
ageAdj_dat[,1] <- NULL
ageAdj_dat <- data.matrix(t(ageAdj_dat))

phen_fi <- "/home/alastairm/Desktop/Landstrom/data/oncollAnalysis_ageAdj_DDS_phen.csv"                          ## file path
phen_dat <- read.table(phen_fi, sep=",", header=TRUE, check.names=FALSE)
rownames(phen_dat) <- phen_dat[,1]                                                                          ## set rownames
phen_dat[,1] <- NULL

## rm 0
ageAdj_dat[ageAdj_dat<0] <- 0
#ageAdj_dat <- apply (ageAdj_dat, c (1, 2), function (x) {
#  (as.integer(x))
#})

d <- DGEList(counts=ageAdj_fi, group=rep(1:2,each=2))
dim(d)
colnames(d)


y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
d <- DGEList(counts=y, group=rep(1:2,each=2))
