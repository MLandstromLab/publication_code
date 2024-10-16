# --- FUNCTIONS FOR PLOTTING KM BY COPY-NUMBER STATUS --- #

# Helper function 
groupByCNstatus = function(cn_status){
  if (cn_status == 0) {
    return("Neutral")
  } else if (cn_status == -1) {
    return("Deleted")
  } else {
    return("Amplified")
  }
} 

#
# Plot KM. Groups based on copy-number status into two groups 
# (Neutral or group of interest defined by mutation param)
#
plotKMcn = function(variable,
                    training_data,
                    end_point_event,
                    end_point_time,
                    mutation,
                    plot_title,
                    life.table = T){
  
  # Subset training data by variables 
  trainingDataFinal =  training_data %>% 
    dplyr::select(all_of(c(variable, end_point_event,  end_point_time)))
  
  # Divide into groups based on CN status 
  
  # Divide into groups 
  group = map_chr(trainingDataFinal[[variable]], groupByCNstatus)
  
  # Add as variable 
  trainingDataFinal$group = group
  
  # Select individuals with Neutral or mutation of interest 
  trainingDataFinal = dplyr::filter(trainingDataFinal, group %in% c("Neutral", mutation))
  
  # 
  if (sum(trainingDataFinal$group == mutation) != 0) {
    
    # Generate the formula for the model 
    survExpression = paste0("Surv(", end_point_time, ", " , end_point_event, ")")
    f <- as.formula(paste(survExpression, "group", sep = " ~ "))
    sFit <- surv_fit(f, data =  trainingDataFinal)
    
    # Store plot and p-value 
    sFit.res = list()
    
    # Prepare KM
    if (life.table == T) {

        sFit.res$Plot <- ggsurvplot(sFit, 
                             data = trainingDataFinal, legend = "bottom",
                             title = paste0(plot_title, " (n = ", nrow(trainingDataFinal) ,")"),
                             legend.title = "Cn status group", pval = TRUE, xlab = "Time (days)", 
                             font.family = "Helvetica", font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15,
                             conf.int = T,
                             break.time.by = 250,
                             surv.plot.height = 0.65, risk.table = TRUE, cumevents = F,
                             ggtheme = theme_classic(), 
                             fontsize = 5, pval.size = 7, tables.font.tickslab = 5, tables.y.text.col = T, tables.y.text = FALSE, 
                             tables.theme = theme_survminer())
    } else {

        sFit.res$Plot <- ggsurvplot(sFit, 
                             data = trainingDataFinal, legend = "none",
                             title = paste0(plot_title, " (n = ", nrow(trainingDataFinal) ,")"),
                             legend.title = "Cn status group", pval = TRUE, xlab = "Time (days)", 
                             font.family = "Helvetica", font.x = 8, font.y = 8, font.tickslab = 8, font.legend = 15,
                             conf.int = T,
                             break.time.by = 1000,
                             surv.plot.height = 0.65, risk.table = FALSE, cumevents = F,
                             ggtheme = theme_classic(), fontsize = 5, pval.size = 5)
    }


    sFit.res$Pval = surv_pvalue(sFit)
    
    return(sFit.res)
  } else {
    # Store plot and p-value 
    sFit.res = list()
    sFit.res$Plot = NULL
    sFit.res$Pval = NULL 
    return(sFit.res)
  }
}

# --- FUNCTION FOR PLOTTING KM BY EXPRESSION STATUS --- # 

#
# Main function
#
plotKMexp = function(variable,
                     training_data,
                     end_point_event,
                     end_point_time,
                     plot_title,
                     life.table = T){
  
  # Subset training data by variables 
  trainingDataFinal =  training_data %>% 
    dplyr::select(all_of(c(variable, end_point_event,  end_point_time)))
  
  # Calculate median expression 
  medianExpression = median(trainingDataFinal[[variable]])
  
  # Divide into two groups 
  group = ifelse(trainingDataFinal[[variable]] > medianExpression, "High", "Low")
  
  # Add as variable 
  trainingDataFinal$group = group
  
  # Generate the formula for the model 
  survExpression = paste0("Surv(", end_point_time, ", " , end_point_event, ")")
  f <- as.formula(paste(survExpression, "group", sep = " ~ "))
  sFit <- surv_fit(f, data =  trainingDataFinal)
 

  # Store plot and p-value 
  sFit.res = list()

  # Prepare KM
  if (life.table == T) {
    sFit.res$Plot <- ggsurvplot(sFit, 
                           data = trainingDataFinal, legend = "bottom",
                           title = paste0(plot_title, " (n = ", nrow(trainingDataFinal) ,")"),
                           legend.title = "Expression group", pval = TRUE, xlab = "Time (days)", 
                           font.family = "Helvetica", font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15,
                           conf.int = T,
                           break.time.by = 250,
                           surv.plot.height = 0.65, risk.table = TRUE, cumevents = F,
                           ggtheme = theme_classic(), 
                           fontsize = 5, pval.size = 7, tables.font.tickslab = 5, tables.y.text.col = T, tables.y.text = FALSE, 
                           tables.theme = theme_survminer())
  } else {
    sFit.res$Plot <- ggsurvplot(sFit, 
                           data = trainingDataFinal, legend = "none",
                           title = paste0(plot_title, " (n = ", nrow(trainingDataFinal) ,")"),
                           legend.title = "Expression group", pval = TRUE, xlab = "Time (days)", 
                           font.family = "Helvetica", font.x = 10, font.y = 10, font.tickslab = 10, font.legend = 8,
                           conf.int = T,
                           break.time.by = 1000,
                           surv.plot.height = 0.65, risk.table = FALSE, cumevents = F,
                           ggtheme = theme_classic(), 
                           fontsize = 5, pval.size = 5)
  }

  # Get p-value 
  sFit.res$Pval = surv_pvalue(sFit)
  
  return(sFit.res)
  
}

# --- FUNCTIONS FOR RUNNING THE UNIVARIATE KM ANALYSIS FOR all genes --- #

#
# HELPER FUNCTIONS
#

#
# Detect present data types in the data
#
detectDataTypes = function(variables){
    res = map(c(".cn", ".exp", ".met", ".mut"), 
        .f = function(x, var){any(str_ends(var, x))},
        var = variables)
    names(res) = c(".cn", ".exp", ".met", ".mut")
    return(res[res == T])
} 

extractPvalue = function(x){
    if (is.null(x$Pval$pval) == F) {
        return(x$Pval$pval)
    } else {
        return(NA)
    }
}

# Main function running Univariate KM analysis for 
# a list of variables


runUnivariateKM = function(input.data, 
                variables,
                clinical.endpoint,
                out.dir,
                plots){

    # Detect data types automatically based on suffix 
    # so far we only consider copynumber and expression features 
    # .cn = copynumber 
    # .exp = expression (RNA-seq)
    # .met = methylation 
    # .mut = small variants 

    # Remove features that are not present 
    variables = variables[variables %in% colnames(input.data$train)] 

    data.types.present = detectDataTypes(variables)

    # Split variables by types and run the KM analysis 
    variables.df = data.frame(variables) 
    colnames(variables.df) = "Variable"
    variables.df$Dtype = str_extract(variables.df$Variable, regex("\\.\\w+"))

    # Data frame containing variable name and the type (suffix)
    variables.by.dtype = split(variables.df, f = variables.df$Dtype)
    
    # Remove clinical features 
    variables.by.dtype = variables.by.dtype[names(variables.by.dtype) %in% names(data.types.present)]   

    # Temporary variables storing results 
    results.training = NULL
    results.validation = NULL
    result.training.table = NULL 
    result.validation.table = NULL 

    # Store all the results tables for different datatypes into 
    # a list which will be converted to a single data frame 
    result.table.ls = list() 

    # Run analysis for expression features if present in the data
    if (".exp" %in% names(variables.by.dtype)) {
    
        results.training = mclapply(X = variables.by.dtype$.exp$Variable, 
              FUN = plotKMexp, 
              mc.cores = 1,
              training_data = input.data$train,
              end_point_event = paste0(clinical.endpoint, ".clin"),
              end_point_time = paste0(clinical.endpoint, ".time.clin"),
              plot_title = "Training data")
        names(results.training) = variables.by.dtype$.exp$Variable

        # Save the KM plots 
        if (plots == T) {
            out.exp.train.dir <- file.path(out.dir, "Expression/Training_data")
            if(!exists(out.exp.train.dir)) dir.create(out.exp.train.dir, recursive = T)
        
            # Iterate over features 
            for (i in 1:length(results.training)) {

                # Extract gene name 
                var = unlist(strsplit(variables.by.dtype$.exp$Variable[[i]], "\\."))[1]

                # Extract the plot 
                plt = results.training[[i]]$Plot

                # Save the KM plot 
                pdf(file.path(out.exp.train.dir, paste0(var, "_expression_KM_with_training_data.pdf")), 
                    width = 15, height = 12, onefile = F)
                print(plt)
                dev.off()
            }
        }

        # Store the p-values into a table 
        pvalues = map_dbl(results.training, .f = extractPvalue)
        result.training.table = data.frame(Feature = variables.by.dtype$.exp$Variable,
                                           pvalues.training = pvalues)

        results.validation = mclapply(X = variables.by.dtype$.exp$Variable, 
              FUN = plotKMexp, 
              mc.cores = 1,
              training_data = input.data$validation,
              end_point_event = paste0(clinical.endpoint, ".clin"),
              end_point_time = paste0(clinical.endpoint, ".time.clin"),
              plot_title = "Validation data")
        names(results.validation) = variables.by.dtype$.exp$Variable

        # Save the KM plots 
        if (plots == T) {
            out.exp.valid.dir <- file.path(out.dir, "Expression/Validation_data")
            if(!exists(out.exp.valid.dir)) dir.create(out.exp.valid.dir, recursive = T)
        
            # Iterate over features 
            for (i in 1:length(results.validation)) {

                # Extract gene name 
                var = unlist(strsplit(variables.by.dtype$.exp$Variable[[i]], "\\."))[1]

                # Extract the plot 
                plt = results.validation[[i]]$Plot

                # Save the KM plot 
                pdf(file.path(out.exp.valid.dir, paste0(var, "_expression_KM_with_validation_data.pdf")), 
                    width = 15, height = 12, onefile = F)
                print(plt)
                dev.off()
            }
        }

        # Store the p-values into a table 
        pvalues = map_dbl(results.validation, .f = extractPvalue)
        result.validation.table = data.frame(Feature = variables.by.dtype$.exp$Variable,
                                           pvalues.validation = pvalues)

        # Merge training and validation data results
        results.merged = dplyr::left_join(result.training.table, result.validation.table, by = "Feature")

        # Update the result table list 
        result.table.ls$exp = results.merged
    }

    # Run analysis for expression features if present in the data
    if (".cn" %in% names(variables.by.dtype)) {
   
        #
        # Run the analysis comparing amplifield with neutral
        #

        results.training = mclapply(X = variables.by.dtype$.cn$Variable, 
              FUN = plotKMcn, 
              mc.cores = 1,
              training_data = input.data$train,
              end_point_event = paste0(clinical.endpoint, ".clin"),
              end_point_time = paste0(clinical.endpoint, ".time.clin"),
              mutation = "Amplified",
              plot_title = "Training data")
        names(results.training) = variables.by.dtype$.cn$Variable

        # Save the KM plots 
        if (plots == T) {
            out.cn.train.dir <- file.path(out.dir, "CN/Training_data")
            if(!exists(out.cn.train.dir)) dir.create(out.cn.train.dir, recursive = T)
        
            # Iterate over features 
            for (i in 1:length(results.training)) {

                # Extract gene name 
                var = unlist(strsplit(variables.by.dtype$.cn$Variable[[i]], "\\."))[1]

                # Extract the plot 
                plt = results.training[[i]]$Plot

                # Save the KM plot 
                pdf(file.path(out.cn.train.dir, paste0(var, "_cn_amp_KM_with_training_data.pdf")), 
                    width = 15, height = 12, onefile = F)
                print(plt)
                dev.off()
            }
        }

        # Store the p-values into a table 
        pvalues = map_dbl(results.training, .f = extractPvalue)
        result.training.table = data.frame(Feature = variables.by.dtype$.cn$Variable,
                                           pvalues.training = pvalues)

        # Run plotKMexp for validation data 
        results.validation = mclapply(X = variables.by.dtype$.cn$Variable, 
              FUN = plotKMcn, 
              mc.cores = 1,
              training_data = input.data$validation,
              end_point_event = paste0(clinical.endpoint, ".clin"),
              end_point_time = paste0(clinical.endpoint, ".time.clin"),
              mutation = "Amplified",
              plot_title = "Validation data")
        names(results.validation) = variables.by.dtype$.cn$Variable

        # Save the KM plots 
        if (plots == T) {
            out.cn.valid.dir <- file.path(out.dir, "CN/Validation_data")
            if(!exists(out.cn.valid.dir)) dir.create(out.cn.valid.dir, recursive = T)
        
            # Iterate over features 
            for (i in 1:length(results.validation)) {

                # Extract gene name 
                var = unlist(strsplit(variables.by.dtype$.cn$Variable[[i]], "\\."))[1]

                # Extract the plot 
                plt = results.validation[[i]]$Plot

                # Save the KM plot 
                pdf(file.path(out.cn.valid.dir, paste0(var, "_cn_amp_KM_with_validation_data.pdf")), 
                    width = 15, height = 12, onefile = F)
                print(plt)
                dev.off()
            }
        }

        # Store the p-values into a table 
        pvalues = map_dbl(results.validation, .f = extractPvalue)
        result.validation.table = data.frame(Feature = variables.by.dtype$.cn$Variable,
                                           pvalues.validation = pvalues)
        
        # Merge training and validation data results
        results.merged = dplyr::left_join(result.training.table, result.validation.table, by = "Feature")
        
        # Alter the prefix to indicate that feature relates to amplification
        results.merged$Feature = paste(results.merged$Feature, "amp" , sep = ".")
        
        # Update the result table list 
        result.table.ls$cn.amp = results.merged
                                 

        #
        # Run the analysis comparing deleted with neutral
        #

        # Run plotKMcn for training data 
        results.training = mclapply(X = variables.by.dtype$.cn$Variable, 
              FUN = plotKMcn, 
              mc.cores = 1,
              training_data = input.data$train,
              end_point_event = paste0(clinical.endpoint, ".clin"),
              end_point_time = paste0(clinical.endpoint, ".time.clin"),
              mutation = "Deleted",
              plot_title = "Training data")
        names(results.training) = variables.by.dtype$.cn$Variable

        # Iterate over features 
        if (plots == T) {
            for (i in 1:length(results.training)) {

                # Extract gene name 
                var = unlist(strsplit(variables.by.dtype$.cn$Variable[[i]], "\\."))[1]

                # Extract the plot 
                plt = results.training[[i]]$Plot

                # Save the KM plot 
                pdf(file.path(out.cn.train.dir, paste0(var, "_cn_del_KM_with_training_data.pdf")), 
                    width = 15, height = 12, onefile = F)
                print(plt)
                dev.off()
            }
        }

        # Store the p-values into a table 
        pvalues = map_dbl(results.training, .f = extractPvalue)
        result.training.table = data.frame(Feature = variables.by.dtype$.cn$Variable,
                                           pvalues.training = pvalues)
         
        # Run plotKMexp for validation data 
        results.validation = mclapply(X = variables.by.dtype$.cn$Variable, 
              FUN = plotKMcn, 
              mc.cores = 1,
              training_data = input.data$validation,
              end_point_event = paste0(clinical.endpoint, ".clin"),
              end_point_time = paste0(clinical.endpoint, ".time.clin"),
              mutation = "Deleted",
              plot_title = "Validation data")
        names(results.validation) = variables.by.dtype$.cn$Variable

        # Iterate over features 
        if (plots == T) {
            for (i in 1:length(results.validation)) {

                # Extract gene name 
                var = unlist(strsplit(variables.by.dtype$.cn$Variable[[i]], "\\."))[1]

                # Extract the plot 
                plt = results.validation[[i]]$Plot

                # Save the KM plot 
                pdf(file.path(out.cn.valid.dir, paste0(var, "_cn_del_KM_with_validation_data.pdf")), 
                    width = 15, height = 12, onefile = F)
                print(plt)
                dev.off()
            }
        }

        # Store the p-values into a table 
        pvalues = map_dbl(results.validation, .f = extractPvalue)
        result.validation.table = data.frame(Feature = variables.by.dtype$.cn$Variable,
                                           pvalues.validation = pvalues)
        
        # Update the result table list 
        # Merge training and validation data results
        results.merged = dplyr::left_join(result.training.table, result.validation.table, by = "Feature")
       
        # Alter the prefix to indicate that feature relates to deletion 
        results.merged$Feature = paste(results.merged$Feature, "del" , sep = ".")

        result.table.ls$cn.del = results.merged
    } else {
        print("skip")
    }

    # Add here the procedure for mutation and methylation data when needed 
    #
    #

    # Finalise the p-value tables 
    result.table.final = do.call("rbind", result.table.ls)
    return(result.table.final)
}

## --- FUNCTION FOR EXTRACTING SIGNIFICANT FEATURES --- #
getSignificantFeatures = function(x, pvalue.thresh) {
    
    # Find the significant features 
    x.filtered = x %>% dplyr::filter(pvalues.training < pvalue.thresh)
    
    # Return vector of features. Remove amp and del suffixes 
    x.filtered.features =  unlist(map_if(x.filtered$Feature, 
                                         str_detect(x.filtered$Feature, ".cn."), 
                                         function(x){str_remove(x, ".amp|.del")}))

    # Remove possible duplicates 
    x.filtered.features = unique(x.filtered.features)

    return(x.filtered.features)

}


# --- FUNCTION FOR PLOTTING KM BY PREDICTED RISK --- #

plotKMbyRelativeRisk = function(data, rel.risk) {

    # Assign relative risk 
    data$pred <- rel.risk
    data$group <- ifelse(data$pred  < median(data$pred), "Low", "High")

    if (length(table(data$group)) == 2) {

        # Fit the survival function 
        s.fit.train <- survfit(Surv(time, status) ~group, data = data)

        # Generate the formula for the model 
        survExpression = paste0("Surv(", "time", ", " , "status", ")")
        f <- as.formula(paste(survExpression, "group", sep = " ~ "))
        sFit <- surv_fit(f, data =  data)

        # Initialise results object 
        sFit.res = list()

    # Prepare plot  
        sFit.res$Plot <- ggsurvplot(s.fit.train, 
                           data = data, legend = "bottom",
                           title = paste0("Risk groups in training cohort", " (n = ", nrow(data) ,")"),
                           legend.title = "Risk group", pval = TRUE, xlab = "Time (days)", 
                           font.family = "Helvetica", font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15,
                           conf.int = T,
                           break.time.by = 250,
                           surv.plot.height = 0.65, risk.table = TRUE, cumevents = F,
                           ggtheme = theme_classic(), 
                           fontsize = 5, pval.size = 7, tables.font.tickslab = 5, tables.y.text.col = T, tables.y.text = FALSE, 
                           tables.theme = theme_survminer())

        # Store relevant information as table 
        sFit.res$table = data.frame(N = nrow(data), 
                                Pvalue = surv_pvalue(sFit))

    } else {
        # Based on the predicted risk it was impossible to group into two groups
        # we cannot get a p-value 
        sFit.res = NULL
    } 
    return(sFit.res)
}
