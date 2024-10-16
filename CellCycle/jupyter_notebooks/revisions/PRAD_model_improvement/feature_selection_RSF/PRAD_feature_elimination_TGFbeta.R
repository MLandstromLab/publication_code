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

library(pec)
library(randomForestSRC)

#
# Function for tuning model parameters
#
#

tune_model <- function(data, formula, ntree, seed) {

  set.seed(seed)
  tuning <- tune(formula, data, mtryStart = 2, ntreeTry = ntree, doBest = TRUE)
  node <- tuning$optimal[1]
  mtry <- tuning$optimal[2]

  error <- tuning$rf$err.rate[500]
  res <- tibble(nodesize = node, mtry = mtry, OOB_error = error)

  return(res)
}

#
# Main function for generating the RSF model
# based on training set
#
generateRSFModel = function(training.data,
                            variables,
                            end.point.event,
                            end.point.time,
                            model.name,
                            ntrees = c(50,100,200,500,1000),
                            seed) {


  # Store relevant data into an object
  rsf.obj = list("model" = model.name)

  # Select only variables of interest
  training.data = training.data  %>% dplyr::select(all_of(c(variables, end.point.event, end.point.time)))

  # Generate the formula for the model
  surv.expression = paste0("Surv(", end.point.time, ", " , end.point.event, ")")
  f <- as.formula(paste(surv.expression, ".", sep = " ~ "))

  # Store formula
  rsf.obj$formula = f

  # First tune the parameters
  tune.res <- list()
  for(n in ntrees) {
    tune.res[[as.character(n)]] <- tune_model(training.data, f, n, seed) # Store optimal values for each ntree
  }

  # Generate a table
  tune.res.df = dplyr::bind_rows(tune.res, .id = "ntree")

  # Store tuning results
  rsf.obj$tuning = tune.res.df

  idx.min <- which.min(sapply(tune.res, function(x) x$OOB_error)) # Get ntree with smallest error
  ntree.optm <- ntrees[idx.min]
  nodesize.optm <- tune.res[[idx.min]]$nodesize # Get the best nodesize for the chosen ntree
  mtry.optm <- tune.res[[idx.min]]$mtry # Get the best mtry for the chosen ntree

  # Finally fit the model with optimised parameters
  set.seed(seed)
  rf.fit.train <- list(eval(bquote (rfsrc(.(f), training.data,
                                          ntree = ntree.optm,
                                          mtry = mtry.optm,
                                          nodesize = nodesize.optm,
                                          importance = TRUE,
                                          block.size = 1,
                                          seed = seed))))

  # Store the fitted model
  rsf.obj$fit.model = rf.fit.train

  return(rsf.obj)
}

#
# Iterative feature elimiation
#

iterativeFeatElim = function(data, 
                             variables, 
                             end.point.event,
                             end.point.time,
                             model.name = "model2",
                             drop.out.perc = 0.2,
                             ntrees = c(50,100,200,500,1000),
                             seed){


  # Store results for testing 
  vimp.result.ls = list() 

  # Errors 
  oob.errors = list()
  
  # Num variables
  num.variables = list()
  
  # Store optimal parameters 
  params.ls = list("ntree" = list(),
                   "mtry" = list(),
                   "nodesize" = list())
  
  # Initialise the number of variables
  variables.left = variables[-c(1,2)]
  
  # Iteration
  i = 1
    
  # If the number of variables is less than 0.2 of the original 
  # data stop iteration
  n.var.stop = round(length(variables) * drop.out.perc)
  
  # Iterate 
  while (length(variables.left) > n.var.stop) {
   
    print(paste0("Iteration : ", i))

    # Create model 
    rsf.model = generateRSFModel(data,
                 variables = variables.left,
                 end.point.event = end.point.event,
                 end.point.time = end.point.time,
                 model.name = "model2",
                 ntrees = ntrees,
                 seed = seed)
    
    # Store the params 
    params.ls[["ntree"]][i] = rsf.model$fit.model[[1]]$ntree
    params.ls[["mtry"]][i] = rsf.model$fit.model[[1]]$mtry
    params.ls[["nodesize"]][i] = rsf.model$fit.model[[1]]$nodesize
  
    # Add the oob.error to list   
    oob.errors[[i]] = rsf.model$fit.model[[1]]$err.rate[rsf.model$fit.model[[1]]$ntree]
  
    # Add the number of variables 
    num.variables[[i]] = length(variables.left)
  
    # Calculate VIMP
    vimp.res = vimp(rsf.model$fit.model[[1]])
    
    # Prepare a table
    vimp.table = as.data.frame(vimp.res$importance) %>% 
                      tibble::rownames_to_column("Variable") 
    colnames(vimp.table) = c("Variable", "Importance")
    
    # Order by decreasing Importance such that the worst features
    # go bottom of the table
    vimp.table = dplyr::arrange(vimp.table, desc(Importance))

    # Store results 
    vimp.result.ls[[i]] = vimp.table
    
    # Remove the worst performing features 
    variables.left = vimp.table$Variable[1:((1 - drop.out.perc) * length(variables.left))]
    
    # Increment 
    i = i + 1
  }
  
  # Combine to table
  final.df = data.frame("ntree" = unlist(params.ls$ntree),
             "mtry" = unlist(params.ls$mtry),
             "nodesize" = unlist(params.ls$nodesize),
             "num_variable" = unlist(num.variables),
             "error" = unlist(oob.errors))


  
  # Find row with minumum error 
  min.idx = which(final.df$error == min(final.df$error))
  print(min.idx)
  
  # Select the last row having error equal to minimum. We want to minimize the 
  # number of features
  if (length(min.idx) > 1) {
    min.idx = max(min.idx)
  }
  print(min.idx)
  
  # Optimal parameters 
  optimal.set.param = final.df[min.idx,]
  optimal.set.variables = vimp.result.ls[[min.idx]]
  
  
  # Final results
  final.results = list("Errors" = final.df, 
                       "Features" = vimp.result.ls,
                       "Optimal.param" = optimal.set.param,
                       "Optimal.var" = optimal.set.variables)
  
  return(final.results)

}


#
# MAIN 
#

# Path to rdata
rdata.path = "/lustre/projects/landstrom_core/data/rdata/manuscript_work/revisions/tgfbeta/PRAD/"

# Read in preprocessed data 
tcga.dataset.merged = readRDS(file.path(rdata.path, "tcga.dataset_merged.rds"))

variables.selected.filtered.ls = readRDS("/lustre/projects/landstrom_core/results/prognostic_model_development_revised_cell_cycle/tgfbeta/PRAD/significant.features.ls")

feature.elim.pfi = iterativeFeatElim(data = tcga.dataset.merged$PFI$train,
                  variables = variables.selected.filtered.ls[["PFI"]],
                  end.point.event = "PFI.clin",
                  end.point.time = "PFI.time.clin",
                  model.name = "model2",
                  drop.out.perc = 0.2,
                  ntrees = c(50,100,200,500,1000),
                  seed = 42)  

saveRDS(feature.elim.pfi, file.path(rdata.path, "Feature_elimination_res_PFI.rds"))
