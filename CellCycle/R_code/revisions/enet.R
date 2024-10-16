#
# Function for finding the optimal lambda based on C-index
#

#findOptimalLambda = function(x, y, out.dir, seed = 42){

    # Set the seed 
#    set.seed(seed)

    # 10-fold cross validation 
#    cv.fit <- cv.glmnet(x, 
#                        y, 
#                        family = "cox", 
#                        type.measure = "C", 
#                        maxit = 100000)

    # Store the values into a list 
#    results = list(model = cv.fit,
#                   lambda.min = cv.fit$lambda.min,
#                   lambda.1se = cv.fit$lambda.1se)

#    # Returns the coefficients of the solution corresponding to lambda 
#    est.coef = coef(cv.fit, s = cv.fit$lambda.min) 
#    est.coef <- as.data.frame(as.matrix(est.coef))
#    colnames(est.coef) <- "coef"
#    # Convert coeffiecients to hazard ratios 
#    est.coef$HR <- exp(est.coef$coef)
#    active.k = which(est.coef$coef != 0)
#    active.k.vals = est.coef[active.k,] 
#    # Store the active coefficients
#    results$active.k.vals = active.k.vals 
#
#    # Produce some visualisations 
#    png(file.path(out.dir, "CV10.png"))
#    plot(cv.fit)
#    dev.off()
#
#    return(results)
#}


findOptimalLambda = function(x, y, out.dir, seed = 42){

    # Set the seed 
    set.seed(seed)

    # 10-fold cross validation 
    out = tryCatch(
    {
        cv.fit <- cv.glmnet(x, 
                        y, 
                        family = "cox", 
                        type.measure = "C", 
                        maxit = 100000)

        # Store the values into a list 
        results = list(model = cv.fit,
                   lambda.min = cv.fit$lambda.min,
                   lambda.1se = cv.fit$lambda.1se)

        # Returns the coefficients of the solution corresponding to lambda 
        est.coef = coef(cv.fit, s = cv.fit$lambda.min) 
        est.coef <- as.data.frame(as.matrix(est.coef))
        colnames(est.coef) <- "coef"
        # Convert coeffiecients to hazard ratios 
        est.coef$HR <- exp(est.coef$coef)
        active.k = which(est.coef$coef != 0)
        active.k.vals = est.coef[active.k,] 
        # Store the active coefficients
        results$active.k.vals = active.k.vals 
        return(results)
    },
    error=function(cond) {
        message("Something went wrond during optimisation")
        return(NULL)
    }
    )
    return(out)
}

#
# Function for returning 
#

findCoefForLambda = function(x, y, out.dir, seed = 42){

    # Set the seed 
    set.seed(seed)

    # 10-fold cross validation 
    cv.fit <- cv.glmnet(x.train.mat, 
                        y.train, 
                        family = "cox", 
                        type.measure = "C", 
                        maxit = 100000)

    # Store the values into a list 
    results = list()
                  
    # Return coef lists for all variables 
    coefs.df = as.data.frame(as.matrix(coef(cv.fit, cv.fit$lambda)))
    active.coef.df = coefs.df != 0 

    colnames(coefs.df) = cv.fit$lambda
    colnames(active.coef.df) = cv.fit$lambda

    results$coefs = coefs.df 
    results$active.coefs = active.coef.df
    return(results) 
} 


#
# Function for running the KM analysis for all of 
# lambda 
#

runKMforAll = function(pcox.fit, dir.res.pcox){

    # Store the results into a list 
    results = list(tables = list(), plots = list())
    for (i in 1:length(pcox.fit$model$lambda)) {

        # Predictions for the validation data
        pred.valid <- predict(pcox.fit$model, 
                      newx = x.valid.mat, 
                      s = pcox.fit$model$lambda[i], 
                      type = "response")

        # Fitted relative risk
        rel.risk.valid <- pred.valid[,1] 

        # Stratify validation data into two groups based on the fitted relative risk
        y.data.valid <- as.data.frame(as.matrix(y.valid))

        # Plot KM and extract the p-value  
        KM.valid.by.risk = plotKMbyRelativeRisk(data = y.data.valid, 
                                        rel.risk = rel.risk.valid )

       results$tables[[i]] = KM.valid.by.risk$table 
       results$plots[[i]] = KM.valid.by.risk$Plot 
        
    }
    return(results)
}

#
# MAIN function 
#

runModelDiagnostics = function(pcox.fit, dir.res){

    # Run KM for all lambdas 
    results.all.km.pred = runKMforAll(pcox.fit, dir.res)
    results.all.km.pred.df = do.call("rbind", results.all.km.pred$tables)

    # Add the actual lambda value 
    results.all.km.pred.df$lambda = pcox.fit$model$lambda 

    # Add the mean C-index from the cross-validation  
    results.all.km.pred.df$cvm = pcox.fit$model$cvm

    # Write to a table 
    write.csv(results.all.km.pred.df, file.path(dir.res, "Diagnostics_all_lambda_KM_evaluation.csv"))

    #We will first explore the relationship betwen the observed P-value from the KM validation and 
    #the mean C-index for each lambda.
    gg = ggplot(results.all.km.pred.df, aes(Pvalue.pval, cvm)) + 
                geom_point() + theme_bw()
    ggsave(plot = gg, filename = file.path(dir.res, "Diagnostics_pvalue_vs_training_cindex.pdf"))

    gg = ggplot(results.all.km.pred.df, aes(lambda, Pvalue.pval)) + 
        geom_line() + theme_bw()
    ggsave(plot = gg, filename = file.path(dir.res, "Diagnostics_lambda_vs_pvalue.pdf"))
}
