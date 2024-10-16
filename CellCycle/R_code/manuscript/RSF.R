library(pec)
library(randomForestSRC)
library(survival)

# splitCases --------------------------------------------------------------

#
# Function for splitting data randomly into 
# training and validation set 
#
splitCases = function(data, split, only.complete = T, variables, seed, impute.valid = F, impute.var = NULL) {
  
  # Number of cases in total
  ncases.initial = nrow(data)
  
  if (only.complete == T) {
     print(variables)
     print(data[,variables])
     # Take only complete cases if data contains missing values 
     print("Taking only complete cases")
     complete <- complete.cases(data[,variables]) # Indexes of cases with no missing values among input variables
     data <- data[complete, variables]
     print(paste0("Including ", sum(complete)  ," cases out of ", ncases.initial, " cases"))
  } else {
     # Allow missing values  	  
     print("Taking all cases")
     data <- data[,variables]
  }
  
  
  if (split != 1.0) {
    
    set.seed(seed)
    trainIdx <- sample(nrow(data), split*nrow(data), replace = FALSE)
    data.train <- data[trainIdx,]
    data.valid <- data[-trainIdx,]
    if (impute.valid == T) {
	print("Performing imputation for validation data")
	set.seed(seed)

    	# Select only features of interest for imputation
	data.imp <- randomForestSRC::impute(data = dplyr::select(data.valid, all_of(impute.var)))
	data.oth = dplyr::select(data.valid, !all_of(impute.var)) 
	data.valid = data.frame(data.oth, data.imp, check.names = F)
    }
    
  } else {
    
    data.train <- data
    data.valid <- data
    
  } 
  result.split = list("train" = data.train,
                      "validation" = data.valid)
  return(result.split)
}


# generateRSFMode ---------------------------------------------------------

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
                                          seed = -13))))
  
  # Store the fitted model
  rsf.obj$fit.model = rf.fit.train
  
  return(rsf.obj)
}


# variableImportance ------------------------------------------------------

variableImportance = function(rsf.obj,
                              dir.results,
                              seed){
  
  # Calculate vimp
  vimp.res = vimp(rsf.obj$fit.model[[1]])
  
  vimp.table = as.data.frame(vimp.res$importance) %>% 
                      tibble::rownames_to_column("Variable") 
  
  colnames(vimp.table) = c("Variable", "Importance")
  vimp.table = dplyr::arrange(vimp.table, desc(Importance))
  
  # Store the table 
  write.csv(vimp.table, file = file.path(dir.results, paste0(rsf.obj$model,"_variable_importance_table.csv")))
  
  # Prepare a plot
  pdf(file.path(dir.results, paste0(rsf.obj$model, "_variable_importance.pdf")))
  plot(vimp.res)
  dev.off()
  
  
}

# evaluateModel -----------------------------------------------------------

#
# Function for evaluating the model. Calculates C-Index and 
# also produces plots 
#
evaluateModel = function(rsf.obj, 
                         validation.data, 
                         dir.results,
                         seed){
  
  # Store results to an object 
  eval.results = list()
  
  
  # Find the prediction error
  predError <- pec(rsf.obj$fit.model[[1]],
                   formula = rsf.obj$formula, 
                   cens.model = "marginal",
                   data = validation.data)
  
  # Store 
  eval.results$pred.err = predError
  
  # Calculate concordance indices of the models for the validation set
  set.seed(seed)
  Cindex.obj <-  pec::cindex(rsf.obj$fit.model[[1]], 
                             formula = rsf.obj$formula, 
                             cens.model = "marginal", 
                             data = validation.data, 
                             confInt = TRUE)
  
  # Store 
  eval.results$Cindex = Cindex.obj$AppCindex[[1]]
  return(eval.results)
}


# iterativeFeatElim -----------------------------------------------------------
# Function performs iterative feature elimination. At each iteration worst features 
# are removed. The percentage of features to remove is given by the user. At each iteration
# the parameters are tuned.
#

# MAIN function
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
  variables.left = variables
  
  # Iteration
  i = 1
  
  # Iterate 
  while (length(variables.left) > 10) {
   
    print(paste0("Iteration : ", i))	

    # Create model 
    rsf.model = generateRSFModel(data,
                 variables = variables.left,
                 end.point.event = "PFI.clin",
                 end.point.time = "PFI.time.clin",
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



