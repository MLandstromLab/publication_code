setwd("/home/data/project_code/landstrom_core/prognostic_model_development_revision/r/notebooks")

library(tidyverse)
library(survival)
library(survminer)
library(glmnet)
library(WriteXLS)
library(ggfortify)
library(circlize)
library(ComplexHeatmap)
library(parallel)
library(broom)
library(survcomp)
library(survivalROC)
library(gtsummary)
source("../getTCGAData.R")
source("../preprocessTCGAData.R")
source("../KM_analysis.R")
source("../Heatmaps.R")
source("../enet.R")

# Number of de genes 
top.de = 150

# Define the cancer type 
cancer.type = "PRAD"

# Read in the table including the clinical features for each cancer type
clin.feat.tb = read.table("/lustre/projects/landstrom_core/data/clin_features_final.csv", sep = "\t", header = T)

# Get Clinical variables
clin.var = unlist(strsplit(clin.feat.tb$Features[clin.feat.tb$Ctype == cancer.type], split = ","))

# Ensembl id mapping file 
ens.id.mapping = "/home/organisms/Human/hg38/Homo_sapiens.GRCh38_March2022/ENSEMBLE_to_SYMBOL.csv"

# Output dir 
out.dir.data = file.path("/lustre/projects/landstrom_core/data/rdata/manuscript_work/revisions", cancer.type)
dir.create(out.dir.data, recursive = T)

clin.var

tcga.cn = getTCGACopyNumberData(cancer.type)

tcga.expr = getTCGAExpressionData(cancer.type, annotation.file = ens.id.mapping)

# Get cancer specific clinical data 
tcga.clin = getClinData(cancer.type)

# Get the end point related clinical data 
tcga.endpoints = getClinEndpointData(cancer.type) %>% dplyr::select(bcr_patient_barcode, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time)

# Merge end point data to clinical data 
tcga.clin = dplyr::left_join(tcga.clin, tcga.endpoints, by = "bcr_patient_barcode")

write.csv(tcga.clin, file.path(out.dir.data, "clinical_data.csv"))

tcga.cn = selectPrimaryT(tcga.cn)

tcga.cn = dropDuplicateSamples(tcga.cn)

tcga.cn.datamat = prepareDataMatrix(tcga.cn)
saveRDS(tcga.cn.datamat, file = file.path(out.dir.data, "copy_number_status.rds"))

tcga.expr = selectPrimaryT(tcga.expr)

tcga.expr = dropDuplicateSamples(tcga.expr)

tcga.expr.datamat = prepareDataMatrix(tcga.expr)

saveRDS(tcga.expr.datamat, 
        file = file.path(out.dir.data, "raw_expressions.rds"))

tcga.dataset = mergeTCGAdata(clin.data = tcga.clin,
                                  data = list("CN" = tcga.cn.datamat, 
                                              "EXP" = tcga.expr.datamat), 
                                  data.suffixes = c("cn","exp"))

# Define function for adding the clinical variables 
addClinVar = function(data, clin.var) {
    if ("Age" %in% clin.var) {
        data$CLIN$Age.clin <- data$CLIN$age_at_diagnosis.clin
    } 
    if ("Tumor.stage" %in% clin.var){
        data$CLIN$Tumor.stage.clin = factor(map_chr(data$CLIN$ajcc_pathologic_stage.clin, reformatTumorStage))
    }
    if ("Gender" %in% clin.var){
        data$CLIN$Gender.clin <- factor(data$CLIN$gender.clin)    
    } 
    if ("Gleason.group" %in% clin.var) {
        
        # Determine the Gleason group 
        data$CLIN$Gleason.group.clin = map2_chr(data$CLIN$primary_gleason_grade.clin, 
                                           data$CLIN$secondary_gleason_grade.clin, 
                                           determineGleasonGroup)

        # Set up the factor levels 
        data$CLIN$Gleason.group.clin = factor(data$CLIN$Gleason.group.clin, 
                                    levels = c("Gleason_group_1", "Gleason_group_2"))
    }
    return(data)
}

# Add clinical variables to dataset
tcga.dataset = addClinVar(tcga.dataset, clin.var)

saveRDS(tcga.dataset, file = file.path(out.dir.data, "tcga.dataset.rds"))

rm("tcga.expr")
rm("tcga.cn")
rm("tcga.clin")

# Read in the preprocessed dataset if continued 
#tcga.dataset = readRDS(file.path(out.dir.data, "tcga.dataset.rds"))

# Raw expression data 
#tcga.expr.raw.datamat = readRDS(file.path(out.dir.data, "raw_expressions.rds"))

# Define and create the root directory for results 
dir.res.root = file.path("/lustre/projects/landstrom_core/results/prognostic_model_development_revised_cell_cycle/DE_based/top150", cancer.type)
dir.create(dir.res.root, recursive = T)

# Define and create the results for the KM analysis 
dir.res.km = file.path(dir.res.root, "Kaplan_Meier_plots")
dir.create(dir.res.km)

# Load the DE table 
de.pfi = read.csv("/lustre/projects/landstrom_core/results/PRAD_specific_analysis_revised/DE_by_predicted_risk/PRAD/PFI_by_risk_de_results_lfc_05.csv")
de.os = de.pfi
#de.os = read.csv("/lustre/projects/landstrom_core/results/PRAD_specific_analysis_revised/DE_by_predicted_risk/PRAD/OS_by_risk_de_results_lfc_05.csv")
de.dfi = read.csv("/lustre/projects/landstrom_core/results/PRAD_specific_analysis_revised/DE_by_predicted_risk/PRAD/DFI_by_risk_de_results_lfc_05.csv")
de.dss = read.csv("/lustre/projects/landstrom_core/results/PRAD_specific_analysis_revised/DE_by_predicted_risk/PRAD/DSS_by_risk_de_results_lfc_05.csv")

# Read in the original gene list 
customer.genes.df = read.csv("/lustre/projects/landstrom_core/data/Customer_genes.tsv", sep = "\t")
customer.genes = customer.genes.df$V1

de.pfi = unique(c(de.pfi$X[1:top.de], customer.genes))
de.os = unique(c(de.os$X[1:top.de], customer.genes))
de.dfi = unique(c(de.dfi$X[1:top.de], customer.genes))
de.dss = unique(c(de.dss$X[1:top.de], customer.genes))


# Test with 150 
gene.lists = list(PFI = de.pfi, 
                  DFI = de.os, 
                  OS = de.dfi,
                  DSS = de.dss)


clinical.end.point.stats = tcga.dataset$CLIN %>% 
                                   dplyr::select(c("OS.clin","DSS.clin","DFI.clin","PFI.clin")) %>%
                                   pivot_longer(everything()) %>%
                                   mutate(value = factor(value)) %>%
                                   group_by(name, value) %>%
                                   summarise(N = n()) %>% 
                                   pivot_wider(names_from =  value,
                                               values_from = N)

# Here we store all the training and validation splits 
train_and_validation_ls = list()

# Variables selected 
variables_selected_ls = list()

# Number of samples in training and validation cohorts 
nsamples_step1_ls = list()

selectVariables = function(clinical.endpoint,
                           clinical.variables,
                           gene.list, 
                           data.suffixes){

    # Construct the clinical features   
    clinical.features =  c(paste0(clinical.endpoint, ".clin"),
                           paste0(clinical.endpoint, ".time.clin"),
                           paste0(clinical.variables, ".clin"))

    # Constructed the sequencing data features 
    seq.features = unlist(map(data.suffixes, 
                           .f = function(x, gene.list){paste(gene.list, x, sep = ".")}, 
                           gene.list = gene.list))

    # Construct the list of variables 
    feature.ls = c(clinical.features, seq.features)

    return(feature.ls)
}



#
# Function for splitting data randomly into 
# training and validation set 
#
splitCases = function(obj, split, variables, seed) {
    
    
  # keep track of statistics 
  ncases.initial = c()
    
  complete.ls = list()
    
  # We will only select the cases with complete data 
  for (i in 1:length(names(obj))){
      
      # Extract the data 
      data.type = names(obj)[i]
      data = obj[[data.type]]
      
      # Find the variables of interest 
      variables.dat = variables[variables %in% colnames(data)]
      
      
      # Check how many individuals there are with non missing data
      # Boolean
      complete = complete.cases(data[,variables.dat])
      
      complete.ls[[i]] = as.data.frame(complete)
      
      # Update the number of complete cases
      ncases.initial = c(ncases.initial, sum(complete))
      
      
  }
    
  # Select complete cases 
  complete.across.all = apply(bind_cols(complete.ls),1, all)
    
  print(paste0("Including ", sum(complete.across.all)  ," cases out of ", max(ncases.initial), " cases"))

  
  # Get the list of complete samples to be included 
  samples.included = obj$CLIN$Participant.ID[complete.across.all]
  
  if (split != 1.0) {
    
    set.seed(seed)
    trainIdx <- sample(length(samples.included), split*length(samples.included), replace = FALSE)
    
    samples.train <- samples.included[trainIdx]
    samples.valid <- samples.included[-trainIdx]
       
  } else {
    
    data.train <- samples.included
    data.valid <- samples.included
    
  } 
  result.split = list("train" = samples.train,
                      "validation" = samples.valid)
  return(result.split)
}



for (end.point in c("OS","DSS","DFI","PFI")){

    
    #
    # No need to change the first part 
    #

    #
    # MODIFICATION. HERE WE SELECT CLIN ENDPOINT SPECIFIC VARIABLES (DE GENES)
    # 
    gene.list = gene.lists[[end.point]]
    
    # Selected variables 
    variables.selected = selectVariables(clinical.endpoint = end.point,
                                         clinical.variables = clin.var,
                                         gene.list = gene.list, 
                                         data.suffixes = c("cn","exp"))
    
    variables_selected_ls[[end.point]] = variables.selected
    
    
    #
    # Splitting function needs to be change to accommondate the 
    # altered data structure 
    #
    
    # Data set is split randomly into training and validation sets. Only complete cases 
    # are selected.
    train_and_validation = splitCases(obj = tcga.dataset, 
                                  split = 0.75, 
                                  variables = variables.selected, 
                                  seed = 42)
    
    # Update list
    train_and_validation_ls[[end.point]] = train_and_validation 
    
    
    # Store number of  
    nsamples.step1 = c(length(train_and_validation$train), length(train_and_validation$validation))
    names(nsamples.step1) = c("ntrain.step1", "nvalid.step1")
    nsamples_step1_ls[[end.point]] = nsamples.step1
}

splitDataset = function(obj, train_and_validation_ls){

    # 
    split.obs.by.endpoint = list()
    
    for (end.point in names(train_and_validation_ls)){
        
        train.samples = train_and_validation_ls[[end.point]]$train
        validation.samples = train_and_validation_ls[[end.point]]$validation
        
        # New entry 
        split.obs.by.endpoint[[end.point]] = list()
        
        for (data.type in names(obj)){
            
            split.obs.by.endpoint[[end.point]][["CLIN"]] = list()  
            split.obs.by.endpoint[[end.point]][["CLIN"]][["train"]] = obj$CLIN %>% 
                          dplyr::filter(Participant.ID %in% train.samples)
            split.obs.by.endpoint[[end.point]][["CLIN"]][["validation"]] = obj$CLIN %>% 
                          dplyr::filter(Participant.ID %in% validation.samples)
                    
            split.obs.by.endpoint[[end.point]][["EXP"]] = list()    
            split.obs.by.endpoint[[end.point]][["EXP"]][["train"]] = obj$EXP %>% 
                          dplyr::filter(Participant.ID %in% train.samples)
            split.obs.by.endpoint[[end.point]][["EXP"]][["validation"]] = obj$EXP %>% 
                          dplyr::filter(Participant.ID %in% validation.samples)
            
            split.obs.by.endpoint[[end.point]][["CN"]] = list()
            split.obs.by.endpoint[[end.point]][["CN"]][["train"]] = obj$CN %>% 
                          dplyr::filter(Participant.ID %in% train.samples)           
            split.obs.by.endpoint[[end.point]][["CN"]][["validation"]] = obj$CN %>% 
                          dplyr::filter(Participant.ID %in% validation.samples)
        }
    }
    return(split.obs.by.endpoint)
}

tcga.dataset.splitted = splitDataset(tcga.dataset, train_and_validation_ls)

convertAge = function(x){
    return(x/360)
}

prepareSummary = function(end.point, data){
    # Convert age 
    data[[end.point]]$CLIN$train$Age.clin = convertAge(data[[end.point]]$CLIN$train$Age.clin)
    data[[end.point]]$CLIN$validation$Age.clin = convertAge(data[[end.point]]$CLIN$validation$Age.clin)
    
    a = data[[end.point]]$CLIN$train %>%
          tbl_summary(include = paste0(clin.var, ".clin")) %>% as_tibble()
    
    b = data[[end.point]]$CLIN$validation %>%
          tbl_summary(include = paste0(clin.var, ".clin")) %>% as_tibble()
    
    test = cbind(a,b)
    test = test[,-3]
    test = t(test)
    test = as.data.frame(test)
    colnames(test) = test[1,]
    test = test[-1,]
    test$End.point = end.point
    return(test)
}

clin.summary.table = bind_rows(map(c("OS","DSS","PFI","DFI"), 
                prepareSummary, 
                data = tcga.dataset.splitted))



# Create dir for debugging 
dir.create("/lustre/projects/landstrom_core/results/debug")

for (end.point in names(train_and_validation_ls)){
    
    # Debug 
    debug.out = file.path("/lustre/projects/landstrom_core/results/debug", end.point)
    dir.create(debug.out)

    # Counts and VST for training data
    counts.training = expDataToMatrix(tcga.dataset.splitted[[end.point]]$EXP$train)
    saveRDS(counts.training, file.path(debug.out, "counts_training.rds"))
    
    vst.transf.training.obj = performVSTtraining(counts.training)
    
    # Counts for evaluation 
    counts.validation = expDataToMatrix(tcga.dataset.splitted[[end.point]]$EXP$validation)
    saveRDS(counts.validation, file.path(debug.out, "counts_validation.rds"))
    vst.transf.validation.counts = performVSTtest(counts = counts.validation, 
                                                  disp.function.train = vst.transf.training.obj$disp.function)
         
    tcga.dataset.splitted[[end.point]]$EXP$train = MatrixToExpdata(vst.transf.training.obj$vst.counts)
    tcga.dataset.splitted[[end.point]]$EXP$validation = MatrixToExpdata(vst.transf.validation.counts)
    
}

saveRDS(tcga.dataset.splitted, file.path(out.dir.data, "tcga.dataset_splitted.rds"))

rm(tcga.dataset)

tcga.dataset.splitted = readRDS(file.path(out.dir.data, "tcga.dataset_splitted.rds"))

tcga.dataset.merged =  mergeDataTypes(tcga.dataset.splitted )

saveRDS(tcga.dataset.merged, file.path(out.dir.data, "tcga.dataset_merged.rds"))

rm(tcga.dataset.splitted)

tcga.dataset.merged = readRDS(file.path(out.dir.data, "tcga.dataset_merged.rds"))

# Store summary statistics 
summary.stats.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    exp.summary.training = prepSummaryExp(x = tcga.dataset.merged[[end.point]]$train, 
                                          raw.data = tcga.expr.datamat,
                                          variables = variables.selected, type = "exp")

    cn.summary.training = prepSummaryCN(tcga.dataset.merged[[end.point]]$train, 
                                        variables = variables.selected, 
                                        type = "cn")
    
    summary.stats.ls[[end.point]] = list("exp.summary.training" = exp.summary.training,
                                         "cn.summary.training" = cn.summary.training)
}


# Store filtered variables 
variables.selected.filtered.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){

    exp.features.keep = summary.stats.ls[[end.point]]$exp.summary.training %>% 
                          filter(`Median expression` > 20, 
                                 `Fraction of zero expression` < 0.75)

    cn.features.keep = summary.stats.ls[[end.point]]$cn.summary.training %>% 
                          filter(`Maximum fraction of aberrations` > 0.15) 

    # Update the summary tables 
    summary.stats.ls[[end.point]]$exp.summary.training$Selected = ifelse(summary.stats.ls[[end.point]]$exp.summary.training$name %in% exp.features.keep$name, "Yes", "No")
    summary.stats.ls[[end.point]]$cn.summary.training$Selected = ifelse(summary.stats.ls[[end.point]]$cn.summary.training$name %in% cn.features.keep$name, "Yes", "No")

    # Collect the variables into vector 
    variables.selected.filtered.ls[[end.point]] = filterFeatures(variables_selected_ls[[end.point]], exp.features.keep$name, type = "exp")
    variables.selected.filtered.ls[[end.point]]= filterFeatures(variables_selected_ls[[end.point]], cn.features.keep$name, type = "cn")
    
}

saveRDS(variables.selected.filtered.ls,  
        file.path(out.dir.data,"variables_selected_filtered_ls.rds"))

variables.selected.filtered.ls = readRDS(file.path(out.dir.data,"variables_selected_filtered_ls.rds"))

# Store the KM tables 
km.pvalue.table.ls = list()

# Store the significant features 
significant.features.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){

    # Create dir for plots 
    dir.create(file.path(dir.res.km, end.point))
    
    if (nrow(tcga.dataset.merged[[end.point]]$train) > 0){
    
        # Run univariate KM
        km.pvalue.table = runUnivariateKM(input.data = tcga.dataset.merged[[end.point]], 
                                          variables = variables.selected.filtered.ls[[end.point]],
                                          clinical.endpoint = end.point,
                                          out.dir = file.path(dir.res.km, end.point),
                                          plots = T)
    
    
        # Sort the results based on the training p-value and write the results to output
        km.pvalue.table = km.pvalue.table %>% dplyr::arrange(pvalues.training)
        km.pvalue.table$Selected = ifelse(km.pvalue.table$pvalues.training < 0.05, "Yes", "No") 
        write.csv(km.pvalue.table, file.path(dir.res.km, paste0(end.point, "_LogRank_pvalues.csv")))
    
        km.pvalue.table.ls[[end.point]] = km.pvalue.table
    
        # Extract the significant features 
        significant.features = getSignificantFeatures(km.pvalue.table, pvalue.thresh = 0.05)

        # Store 
        significant.features.ls[[end.point]] = significant.features
        
    } else {
        significant.features.ls[[end.point]] = NULL
    }
}

saveRDS(significant.features.ls, file.path(dir.res.root, "significant.features.ls"))

dir.res.pcox = file.path(dir.res.root, "Penalized_Cox_risk_prediction/customer_features/Without_clinical_features")
dir.create(dir.res.pcox, recursive = T)

# Helper function for fixing variable names 
fixVarNames = function(x){
    if (str_detect(x, "Gender")) {
        return("Gender")
    } else if (str_detect(x, "Tumor.stage")){
        return("Tumor.stage")
    } else if (str_detect(x,".cn")){
        return(str_extract(x, "\\w+.cn"))
    } else if (str_detect(x, "Gleason.group")){ 
        return("Gleason.group")
    } else {
        return(x)
    }
}

# Number of samples in training and validation cohorts 
nsamples_step2_ls_no_clin = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (is.null(significant.features.ls[[end.point]]) == F) {
    
    
        # Store the number of samples       
        nsamples.step2 = c(nrow(tcga.dataset.merged[[end.point]]$train), 
                           nrow(tcga.dataset.merged[[end.point]]$validation))
    }    
    else {
    
        # If there are no significant features store NULL
        
        # Store 
        nsamples.step2 = c(NA, NA)
    }
        
    names(nsamples.step2) = c("ntrain.step2", "nvalid.step2")
    nsamples_step2_ls_no_clin[[end.point]] = nsamples.step2
}

# Store significant features 
rcox.res.no.clin.ls = list()

# Store model matrices
model.matrices.ls = list()

# Store the fitted models for prediction 
pcox.fit.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    
    # Construct the clinical end points 
    end_point_event = paste0(end.point, ".clin")
    end_point_time = paste0(end.point, ".time.clin")
    
    # Subset 
    selected.features = c(end_point_event, end_point_time, significant.features.ls[[end.point]])
    
    # Input training data 
    input.training = tcga.dataset.merged[[end.point]]$train %>% dplyr::select(all_of(selected.features))
    input.validation = tcga.dataset.merged[[end.point]]$validation %>% dplyr::select(all_of(selected.features))
    
    # Check the number of features
    # Regulariation cannot be run if there is only one feature
    num.feat = ncol(input.training) - 2
    
    if (is.null(input.training) == F){
        if (num.feat > 1) {
    
            # Genereate model matrix 
            model.matrices = generateModelMatrices(input.training, 
                                                   input.validation, 
                                                   clinical.endpoint = end.point)
        
            model.matrices.ls[[end.point]] = model.matrices
    
            # Create output dir 
            dir.create(file.path(dir.res.pcox, end.point))
    
            # Find optimal lambda (hyperparameter for elastic net)
            pcox.fit = findOptimalLambda(x = model.matrices$x.train.mat, 
                             y = model.matrices$y.train,
                             out.dir = file.path(dir.res.pcox, end.point))
        
            pcox.fit.ls[[end.point]] = pcox.fit
    
            # Write the final features included in the model to a file 
            WriteXLS(pcox.fit$active.k.vals, 
            file.path(dir.res.pcox, end.point ,"Active_covariates_in_lambda.min_model.xlsx"), 
            BoldHeaderRow = T,
            row.names = T)
    
            # Final significant features 
            rcox.res.no.clin = pcox.fit$active.k.vals %>% tibble::rownames_to_column("Feature")
            rcox.res.no.clin.ls[[end.point]] = rcox.res.no.clin  
        } else {
            # If no significant features from earlier steps for the clin. end point then store null
            model.matrices.ls[[end.point]] = NULL
            pcox.fit.ls[[end.point]] = NULL
            rcox.res.no.clin.ls[[end.point]] = NULL
        }

    } else {
        # If no significant features from earlier steps for the clin. end point then store null
        model.matrices.ls[[end.point]] = NULL
        pcox.fit.ls[[end.point]] = NULL
        rcox.res.no.clin.ls[[end.point]] = NULL
    }
}

# Store the result tables
KM.train.by.risk.ls = list()

input.train.test = NULL
y.data.test = NULL
pred.train.test = NULL

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (!is.null(pcox.fit.ls[[end.point]])) {
    
        # Predictions for the training set
        pred.train <- predict(pcox.fit.ls[[end.point]]$model, 
                      newx = model.matrices.ls[[end.point]]$x.train.mat, 
                      s = "lambda.min", 
                      type = "response")

        # Fitted relative risk
        rel.risk <- pred.train[,1] 

        # Stratify validation data into two groups based on the fitted relative risk
        y.data <- as.data.frame(as.matrix(model.matrices.ls[[end.point]]$y.train))
        
        
        # Plot KM and extract the p-value  
        KM.train.by.risk = plotKMbyRelativeRisk(data = y.data, 
                                        rel.risk = rel.risk)
        
        
        if (!is.null(KM.train.by.risk)) {
        
            # Store
            KM.train.by.risk.ls[[end.point]] =  KM.train.by.risk$table
    
            # Store the KM plot
            pdf(file.path(dir.res.pcox, end.point ,"glmnet_K-M_plot_with_training_data.pdf"), 
                width = 15, height = 12, onefile = F)
            print(KM.train.by.risk$Plot)
            dev.off()
    
            # Heatmap preparation         
    
            # Variables to be selected 
            # Because Gender has been changed to a dummy variable its name has been changed
            variables.selected = map_chr(rcox.res.no.clin.ls[[end.point]]$Feature, fixVarNames)
            
            # Fix variable names with -
            variables.selected = variables.selected %>% str_replace("`","") %>% str_replace("`","")

            # Get all input variables
            heatmap.input.train = model.matrices.ls[[end.point]]$x.train %>% dplyr::select(all_of(variables.selected))
                
               
            # Heatmap of training data predictions
            hmap.train <- prepareHeatmap(heatmap.input.train, y.data, pred.train, file.path(dir.res.pcox, end.point), "glmnet_training", row.height = 8)          
            
        
        } else {
            KM.train.by.risk.ls[[end.point]] = NULL
        }
        
    } else {
        KM.train.by.risk.ls[[end.point]] = NULL
    }
}

# Store the result tables
KM.valid.by.risk.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (!is.null(pcox.fit.ls[[end.point]])) {
    
        # Predictions for the training set
        pred.valid <- predict(pcox.fit.ls[[end.point]]$model, 
                      newx = model.matrices.ls[[end.point]]$x.valid.mat, 
                      s = "lambda.min", 
                      type = "response")

        # Fitted relative risk
        rel.risk <- pred.valid[,1] 

        # Stratify validation data into two groups based on the fitted relative risk
        y.data <- as.data.frame(as.matrix(model.matrices.ls[[end.point]]$y.valid))

        # Plot KM and extract the p-value  
        KM.valid.by.risk = plotKMbyRelativeRisk(data = y.data, 
                                        rel.risk = rel.risk)
        
        if (!is.null(KM.train.by.risk)) {
        
            # Store
            KM.valid.by.risk.ls[[end.point]] =  KM.valid.by.risk$table
    
            # Store the KM plot
            pdf(file.path(dir.res.pcox, end.point ,"glmnet_K-M_plot_with_validation_data.pdf"), 
            width = 15, height = 12, onefile = F)
            print(KM.valid.by.risk$Plot)
            dev.off()
    
            # Heatmap preparation 
    
            # Variables to be selected 
            variables.selected = map_chr(rcox.res.no.clin.ls[[end.point]]$Feature, fixVarNames) 
            
            variables.selected = variables.selected %>% str_replace("`","") %>% str_replace("`","")
    
            # Get all input variables
            heatmap.input.valid = model.matrices.ls[[end.point]]$x.valid %>% dplyr::select(all_of(variables.selected))
    
            # Heatmap of training data predictions
            hmap.valid <- prepareHeatmap(heatmap.input.valid, y.data, pred.valid, file.path(dir.res.pcox, end.point), "glmnet_validation", row.height = 8)  
        
        } else {
            KM.valid.by.risk.ls[[end.point]] = NULL
        }
        
    } else {
        KM.valid.by.risk.ls[[end.point]] = NULL
    }
}

KM.by.risk.no.clin.train = bind_rows(KM.train.by.risk.ls, .id = "End point")
KM.by.risk.no.clin.valid = bind_rows(KM.valid.by.risk.ls, .id = "End point")

# Store final resilts
write.csv(KM.by.risk.no.clin.train, file.path(dir.res.pcox, "Final_evaluation_results_training.csv"), row.names = F)
write.csv(KM.by.risk.no.clin.valid, file.path(dir.res.pcox, "Final_evaluation_results_validation.csv"),row.names = F)

dir.res.pcox = file.path(dir.res.root, "Penalized_Cox_risk_prediction/customer_features/With_clinical_features")
dir.create(dir.res.pcox, recursive = T)

# Number of samples in training and validation cohorts 
nsamples_step2_ls_with_clin = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (is.null(significant.features.ls[[end.point]]) == F) {
    
    
        # Store the number of samples       
        nsamples.step2 = c(nrow(tcga.dataset.merged[[end.point]]$train), 
                           nrow(tcga.dataset.merged[[end.point]]$validation))
    }    
    else {
    
        # If there are no significant features store NULL
        
        # Store 
        nsamples.step2 = c(NA, NA)
    }
        
    names(nsamples.step2) = c("ntrain.step2", "nvalid.step2")
    nsamples_step2_ls_with_clin[[end.point]] = nsamples.step2
}

# Store significant features 
rcox.res.with.clin.ls = list()

# Store model matrices
model.matrices.ls = list()

# Store the fitted models for prediction 
pcox.fit.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    # Construct the clinical end points 
    end_point_event = paste0(end.point, ".clin")
    end_point_time = paste0(end.point, ".time.clin")
    
    # Clinical features
    clinical.feat = paste0(clin.var, ".clin")

    
    # Subset 
    selected.features = c(end_point_event, end_point_time, clinical.feat, significant.features.ls[[end.point]])
    
    # Input training data 
    input.training = tcga.dataset.merged[[end.point]]$train %>% dplyr::select(all_of(selected.features))
    input.validation = tcga.dataset.merged[[end.point]]$validation %>% dplyr::select(all_of(selected.features))
    
    # Check the number of features
    # Regulariation cannot be run if there is only one feature
    num.feat = ncol(input.training) - 2
    
    if (is.null(tcga.dataset.merged[[end.point]]$train) == F){
        if (num.feat > 1) {
    
            # Genereate model matrix 
            model.matrices = generateModelMatrices(input.training, 
                                                   input.validation, 
                                                   clinical.endpoint = end.point)          
        
            model.matrices.ls[[end.point]] = model.matrices
    
            # Create output dir 
            dir.create(file.path(dir.res.pcox, end.point))
    
            # Find optimal lambda (hyperparameter for elastic net)
            pcox.fit = findOptimalLambda(x = model.matrices$x.train.mat, 
                             y = model.matrices$y.train,
                             out.dir = file.path(dir.res.pcox, end.point))
        
            pcox.fit.ls[[end.point]] = pcox.fit
            
            if (is.null(pcox.fit) == F){
                # Write the final features included in the model to a file 
                WriteXLS(pcox.fit$active.k.vals, 
                     file.path(dir.res.pcox, end.point ,"Active_covariates_in_lambda.min_model.xlsx"), 
                      BoldHeaderRow = T,
                      row.names = T)           
            
                    
                # Final significant features 
                rcox.res.with.clin = pcox.fit$active.k.vals %>% tibble::rownames_to_column("Feature")
                rcox.res.with.clin.ls[[end.point]] = rcox.res.with.clin  
                
            } else {
                model.matrices.ls[[end.point]] = NULL
                pcox.fit.ls[[end.point]] = NULL
                rcox.res.with.clin.ls[[end.point]] = NULL
            
            }
            
        } else {
            # If no significant features from earlier steps for the clin. end point then store null
            model.matrices.ls[[end.point]] = NULL
            pcox.fit.ls[[end.point]] = NULL
            rcox.res.with.clin.ls[[end.point]] = NULL
        }

    } else {
        # If no significant features from earlier steps for the clin. end point then store null
        model.matrices.ls[[end.point]] = NULL
        pcox.fit.ls[[end.point]] = NULL
        rcox.res.with.clin.ls[[end.point]] = NULL
    }
}

# Store the result tables
KM.train.by.risk.ls = list()
C.index.train.ls = list()
AUC.train.ls = list()

# Helper function for fixing variable names 
fixVarNames = function(x){
    if (str_detect(x, "Gender.clin")) {
        return("Gender")
    } else if (str_detect(x, "Tumor.stage.clin")){
        return("Tumor.stage")
    } else if (str_detect(x,".cn")){
        return(str_extract(x, "\\w+.cn"))
    } else if (str_detect(x, "Gleason.group.clin")){ 
        return("Gleason.group.clin")
    } else {
        return(x)
    }
}

# Calculate AUC
calcAUC = function(pred.time, time, status, risk.score){
    
    # Calculate ROC characteristics
    res.ROC = survivalROC(Stime = time,
                          status = status,
                          marker = risk.score,
                          predict.time = pred.time,
                          method  = "KM")
    return(min(res.ROC$AUC,1))
}

# Set seed 
set.seed(42)

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (!is.null(pcox.fit.ls[[end.point]])) {
    
        # Predictions for the training set
        pred.train <- predict(pcox.fit.ls[[end.point]]$model, 
                      newx = model.matrices.ls[[end.point]]$x.train.mat, 
                      s = "lambda.min", 
                      type = "response")

        # Fitted relative risk
        rel.risk <- pred.train[,1] 
        
        # Calculate the C-index (NEW ADDITION )
        c.index.train = Cindex(pred = pred.train, y = model.matrices.ls[[end.point]]$y.train) 
        C.index.train.ls[[end.point]] = c.index.train 

        # Stratify validation data into two groups based on the fitted relative risk
        y.data <- as.data.frame(as.matrix(model.matrices.ls[[end.point]]$y.train))
        
        
        # TEST new function for calculating the C-index
        cindex.train = concordance.index(rel.risk, 
                                         y.data$time, 
                                         y.data$status,
                                         na.rm = TRUE)
        
        C.index.train.ls[[end.point]] = data.frame("C-index" = round(cindex.train$c.index, 4),
                                          "C-index CI" = paste0("(", round(cindex.train$lower, 4), " - ",  
                                                                round(cindex.train$upper, 4), ")"),
                                        check.names = F)
        

        # Plot KM and extract the p-value  
        KM.train.by.risk = plotKMbyRelativeRisk(data = y.data, 
                                        rel.risk = rel.risk)
        
        
        # Calculate AUCs 
        auc.vect = round(map_dbl(c(365, 3 * 365, 5 * 365), .f = calcAUC, 
                   time = y.data$time, status = y.data$status , risk.score = rel.risk), 4)
        names(auc.vect) = c("AUC (1y)", "AUC (3y)", "AUC (5y)")
        
        AUC.train.ls[[end.point]] = as.data.frame(t(data.frame(auc.vect)))
    
        if (!is.null(KM.train.by.risk)) {
        
            # Store
            KM.train.by.risk.ls[[end.point]] =  KM.train.by.risk$table
    
            # Store the KM plot
            pdf(file.path(dir.res.pcox, end.point ,"glmnet_K-M_plot_with_training_data.pdf"), 
                width = 15, height = 12, onefile = F)
            print(KM.train.by.risk$Plot)
            dev.off()
    
            # Heatmap preparation 
    
            # Variables to be selected 
            # Because Gender has been changed to a dummy variable its name has been changed
            variables.selected = map_chr(rcox.res.with.clin.ls[[end.point]]$Feature, fixVarNames)
            variables.selected = variables.selected %>% str_replace("`","") %>% str_replace("`","")
            
    
            # Get all input variables
            heatmap.input.train = model.matrices.ls[[end.point]]$x.train %>% dplyr::select(all_of(variables.selected))
    
            # Heatmap of training data predictions
            hmap.train <- prepareHeatmap(heatmap.input.train, y.data, pred.train, file.path(dir.res.pcox, end.point), "glmnet_training", row.height = 8)
            #print(hmap.train)
            
        } else {
            KM.train.by.risk.ls[[end.point]] = NULL
            C.index.train.ls[[end.point]] = NULL
            AUC.train.ls[[end.point]] = NULL
            
        }
    } else {
        KM.train.by.risk.ls[[end.point]] = NULL
        C.index.train.ls[[end.point]] = NULL
        AUC.train.ls[[end.point]] = NULL
    }
}

# Store the result tables
KM.valid.by.risk.ls = list()
C.index.valid.ls = list()
AUC.valid.ls = list()

# Set seed
set.seed(42)

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (!is.null(pcox.fit.ls[[end.point]])) {
    
        # Predictions for the validation set
        pred.valid <- predict(pcox.fit.ls[[end.point]]$model, 
                      newx = model.matrices.ls[[end.point]]$x.valid.mat, 
                      s = "lambda.min", 
                      type = "response")

        # Fitted relative risk
        rel.risk <- pred.valid[,1] 

        # Stratify validation data into two groups based on the fitted relative risk
        y.data <- as.data.frame(as.matrix(model.matrices.ls[[end.point]]$y.valid))
        
        
        # TEST new function for calculating the C-index
        cindex.valid = concordance.index(rel.risk, 
                                         y.data$time, 
                                         y.data$status,
                                         na.rm = TRUE)
        
        C.index.valid.ls[[end.point]] = data.frame("C-index" = round(cindex.valid$c.index, 4),
                                          "C-index CI" = paste0("(", round(cindex.valid$lower, 4), " - ",  
                                                                round(cindex.valid$upper, 4), ")"),
                                        check.names = F)
        
        # Plot KM and extract the p-value  
        KM.valid.by.risk = plotKMbyRelativeRisk(data = y.data, 
                                        rel.risk = rel.risk)
        
        # Calculate AUCs 
        auc.vect = round(map_dbl(c(365, 3 * 365, 5 * 365), .f = calcAUC, 
                   time = y.data$time, status = y.data$status , risk.score = rel.risk),4)
        
        names(auc.vect) = c("AUC (1y)", "AUC (3y)", "AUC (5y)")
        
        AUC.valid.ls[[end.point]] = as.data.frame(t(data.frame(auc.vect)))
    
        if (!is.null(KM.train.by.risk)) {
            
            # Store
            KM.valid.by.risk.ls[[end.point]] =  KM.valid.by.risk$table
    
            # Store the KM plot
            pdf(file.path(dir.res.pcox, end.point ,"glmnet_K-M_plot_with_validation_data.pdf"), 
                width = 15, height = 12, onefile = F)
            print(KM.valid.by.risk$Plot)
            dev.off()
    
            # Heatmap preparation 
    
            # Variables to be selected 
            variables.selected = map_chr(rcox.res.with.clin.ls[[end.point]]$Feature, fixVarNames)

            variables.selected = variables.selected %>% str_replace("`","") %>% str_replace("`","")

            # Get all input variables
            heatmap.input.valid = model.matrices.ls[[end.point]]$x.valid %>% dplyr::select(all_of(variables.selected))
    
            # Heatmap of training data predictions
            hmap.valid <- prepareHeatmap(heatmap.input.valid, y.data, pred.valid, file.path(dir.res.pcox, end.point), "glmnet_validation", row.height = 8)
            
        } else {
            KM.valid.by.risk.ls[[end.point]] = NULL
            C.index.valid.ls[[end.point]] = NULL
            AUC.valid.ls[[end.point]] = NULL
        }
    } else {
        KM.valid.by.risk.ls[[end.point]] = NULL
        C.index.valid.ls[[end.point]] = NULL
        AUC.valid.ls[[end.point]] = NULL
    }
}

# Collect the results into a single data frame

# Log-rank test results 
KM.by.risk.with.clin.train = bind_rows(KM.train.by.risk.ls, .id = "End point")
KM.by.risk.with.clin.valid = bind_rows(KM.valid.by.risk.ls, .id = "End point")

# C-indices
C.index.train.df = bind_rows(C.index.train.ls , .id = "End point")
C.index.valid.df = bind_rows(C.index.valid.ls , .id = "End point")

# AUCs 
auc.train.df = bind_rows(AUC.train.ls, .id  = "End point") 
auc.valid.df = bind_rows(AUC.valid.ls, .id  = "End point")

# Store final resilts
write.csv(KM.by.risk.with.clin.train, file.path(dir.res.pcox, "Final_evaluation_results_training.csv"), row.names = F)
write.csv(KM.by.risk.with.clin.valid, file.path(dir.res.pcox, "Final_evaluation_results_validation.csv"),row.names = F)

write.csv(C.index.train.df, file.path(dir.res.pcox, "Final_evaluation_C_index_training.csv"), row.names = F)
write.csv(C.index.valid.df, file.path(dir.res.pcox, "Final_evaluation_C_index_validation.csv"), row.names = F)

write.csv(auc.train.df, file.path(dir.res.pcox, "Final_evaluation_AUC_training.csv"), row.names = F)
write.csv(auc.valid.df, file.path(dir.res.pcox, "Final_evaluation_AUC_validation.csv"), row.names = F)

dir.res.pcox = file.path(dir.res.root, "Penalized_Cox_risk_prediction/customer_features/Only_clinical_features")
dir.create(dir.res.pcox, recursive = T)

# Store COX models
pcox.ref.fit.ls = list()

# 
# Function fits a cox regression model
# 
fitCoxModel = function(data, end.point, features){
    
    # Expand to variable name
    end_point_time = paste0(end.point, ".time.clin")
    end_point_event = paste0(end.point, ".clin")

    # Generate a survival formula object 
    survExpression = paste0("Surv(", end_point_time, ", " , end_point_event, ")")
    f <- as.formula(paste(survExpression, paste(features, collapse = " + "), sep = " ~ "))
    
    model.fit = coxph(f, data = data)
    return(model.fit)
}

# Fit the models
for (end.point in c("OS","DSS","DFI","PFI")){
    
    # Construct the clinical end points 
    end_point_event = paste0(end.point, ".clin")
    end_point_time = paste0(end.point, ".time.clin")
    
    clinical.feat = paste0(clin.var, ".clin")
    
    # Subset 
    selected.features = c(end_point_event, end_point_time, clinical.feat)
    
    # Input training data 
    input.training = tcga.dataset.merged[[end.point]]$train %>% dplyr::select(all_of(selected.features))
    
    if (nrow(input.training) > 1) {
        pcox.ref.fit.ls[[end.point]] = fitCoxModel(input.training, end.point, clinical.feat)
    } else {
        pcox.ref.fit.ls[[end.point]] = NULL
    }
}

# Store the result tables
KM.train.ref.by.risk.ls = list()
C.index.ref.train.ls = list()
AUC.train.ref.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    if (is.null(pcox.ref.fit.ls[[end.point]]) != T) {
    
        # Apply model to predict the risk scores 
        rel.risk = predict(object = pcox.ref.fit.ls[[end.point]], 
               newdata = tcga.dataset.merged[[end.point]]$train[,paste0(clin.var, ".clin"), drop = F], 
               type = "risk")

        # Stratify validation data into two groups based on the fitted relative risk
        y.data <- tcga.dataset.merged[[end.point]]$train[paste0(end.point, c(".clin",".time.clin"))]
        colnames(y.data) = c("status","time")
    
    
        # TEST new function for calculating the C-index
        cindex.ref.train = concordance.index(rel.risk, 
                                         y.data$time, 
                                         y.data$status,
                                         na.rm = TRUE)
        
        C.index.ref.train.ls[[end.point]] = data.frame("C-index" = round(cindex.ref.train$c.index, 4),
                                   "C-index CI" = paste0("(", round(cindex.ref.train$lower, 4), " - ",  
                                                         round(cindex.ref.train$upper, 4), ")"), check.names = F)
    
        # Plot KM and extract the p-value  
        KM.train.ref.by.risk = plotKMbyRelativeRisk(data = y.data, 
                                                                rel.risk = rel.risk)
    
        KM.train.ref.by.risk.ls[[end.point]] =  KM.train.ref.by.risk$table
    
        # Calculate AUCs 
        auc.vect = round(map_dbl(c(365, 3 * 365, 5 * 365), .f = calcAUC, 
                   time = y.data$time, status = y.data$status , risk.score = rel.risk),4)
        
        names(auc.vect) = c("AUC (1y)", "AUC (3y)", "AUC (5y)")
        
        AUC.train.ref.ls[[end.point]] = as.data.frame(t(data.frame(auc.vect)))  
    
    } else {
        
        KM.train.ref.by.risk.ls[[end.point]] = NULL
        C.index.ref.train.ls[[end.point]] = NULL
        AUC.train.ref.ls[[end.point]] = NULL
        
    }    
}

# Store the result tables
KM.valid.ref.by.risk.ls = list()
C.index.ref.valid.ls = list()
AUC.valid.ref.ls = list()

# Iterate over end points 
for (end.point in c("OS","DSS","DFI","PFI")){
    
    
    if (is.null(pcox.ref.fit.ls[[end.point]]) != T) {
    
        # Apply model to predict the risk scores 
        rel.risk = predict(object = pcox.ref.fit.ls[[end.point]], 
               newdata = tcga.dataset.merged[[end.point]]$validation[,paste0(clin.var, ".clin"), drop = F], 
               type = "risk")
    
        # Stratify validation data into two groups based on the fitted relative risk
        y.data <- tcga.dataset.merged[[end.point]]$validation[paste0(end.point, c(".clin",".time.clin"))]
        colnames(y.data) = c("status","time")
    
    
        # TEST new function for calculating the C-index
        cindex.ref.valid = concordance.index(rel.risk, 
                                         y.data$time, 
                                         y.data$status,
                                         na.rm = TRUE)
        
        C.index.ref.valid.ls[[end.point]] = data.frame("C-index" = round(cindex.ref.valid$c.index, 4),
                                                   "C-index CI" = paste0("(", round(cindex.ref.valid$lower, 4), " - ",  
                                                         round(cindex.ref.valid$upper, 4), ")"), check.names = F)
    
        # Plot KM and extract the p-value  
        KM.valid.ref.by.risk = plotKMbyRelativeRisk(data = y.data, 
                                                rel.risk = rel.risk)
    
        KM.valid.ref.by.risk.ls[[end.point]] =  KM.valid.ref.by.risk$table
    
    
        # Calculate AUCs 
        auc.vect = round(map_dbl(c(365, 3 * 365, 5 * 365), .f = calcAUC, 
                   time = y.data$time, status = y.data$status , risk.score = rel.risk),4)
        
        names(auc.vect) = c("AUC (1y)", "AUC (3y)", "AUC (5y)")
        
        AUC.valid.ref.ls[[end.point]] = as.data.frame(t(data.frame(auc.vect)))
        
    } else {
        
        KM.valid.ref.by.risk.ls[[end.point]] = NULL
        C.index.ref.valid.ls[[end.point]] = NULL
        AUC.valid.ref.ls[[end.point]] = NULL
    
    }
}

# Collect the results into a single data frame

# Log-rank test results 
KM.train.ref.by.risk.df = bind_rows(KM.train.ref.by.risk.ls, .id = "End point")
KM.valid.ref.by.risk.df = bind_rows(KM.valid.ref.by.risk.ls, .id = "End point")

# C-indices
C.index.ref.train.df = bind_rows(C.index.ref.train.ls , .id = "End point")
C.index.ref.valid.df = bind_rows(C.index.ref.valid.ls , .id = "End point")

# AUCs 
auc.train.ref.df = bind_rows(AUC.train.ref.ls, .id  = "End point") 
auc.valid.ref.df = bind_rows(AUC.valid.ref.ls, .id  = "End point")

# Store final resilts
write.csv(KM.train.ref.by.risk.df, file.path(dir.res.pcox, "Final_evaluation_results_training.csv"), row.names = F)
write.csv(KM.valid.ref.by.risk.df, file.path(dir.res.pcox, "Final_evaluation_results_validation.csv"),row.names = F)

write.csv(C.index.ref.train.df, file.path(dir.res.pcox, "Final_evaluation_C_index_training.csv"), row.names = F)
write.csv(C.index.ref.valid.df, file.path(dir.res.pcox, "Final_evaluation_C_index_validation.csv"), row.names = F)

write.csv(auc.train.ref.df, file.path(dir.res.pcox, "Final_evaluation_AUC_training.csv"), row.names = F)
write.csv(auc.valid.ref.df, file.path(dir.res.pcox, "Final_evaluation_AUC_validation.csv"), row.names = F)

# Collect into a list 
final.result.collection = list("FSS1_sample_summary" = nsamples_step1_ls,
                               "Clinical_endpoint_stats" = clinical.end.point.stats,
                               "Clin_Feature_summary_stats" = clin.summary.table,
                               "Feature_summary_stats" = summary.stats.ls,
                               "Clin_Feature_summary_stats" = clin.summary.table,
                               "FSS1_results_summary" = km.pvalue.table.ls,
                               "FSS2_sample_summary_no_clin" = nsamples_step2_ls_no_clin,
                               "FSS2_regcox_res_no_clin" = rcox.res.no.clin.ls,
                               "Final_res_KM_no_clin_train" = KM.by.risk.no.clin.train,
                               "Final_res_KM_no_clin_valid" = KM.by.risk.no.clin.valid,
                               "FSS2_sample_summary_with_clin" = nsamples_step2_ls_with_clin,
                               "FSS2_regcox_res_with_clin" = rcox.res.with.clin.ls,
                               "Final_res_KM_with_clin_train" = KM.by.risk.with.clin.train,
                               "Final_res_KM_with_clin_valid" = KM.by.risk.with.clin.valid,
                               "FSS2_sample_summary_only_clin" = nsamples.step2,
                               "Final_res_KM_only_clin_train" = KM.train.ref.by.risk.df,
                               "Final_res_KM_only_clin_valid" = KM.valid.ref.by.risk.df)

saveRDS(final.result.collection, file.path(dir.res.root, "Final_results_collection.rds"))

final.result.collection$rcox.res.with.clin.ls


