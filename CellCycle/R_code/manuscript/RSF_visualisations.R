library(survminer)
library(riskRegression)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# -------------------------------------------------------------------------

#
# Diagnostic plots. Typically produced when  
# model is trained
#

diagnosticPlots = function(rsf.obj, 
                           plot.pp = 2,
                           dir.results){
  
  pdf(file.path(dir.results, paste0(rsf.obj$model, "_variables_diagnostic.pdf")))
  plot.variable(rsf.model$fit.model[[1]], 
                plots.per.page = plot.pp)
  dev.off()
  
  pdf(file.path(dir.results, paste0(rsf.obj$model, "_survival_diagnostic.pdf")))
  plot.survival(rsf.model$fit.model[[1]], 
                plots.per.page =  plot.pp)
  dev.off()
  
}

# -------------------------------------------------------------------------

#
# Error plots 
#
#

errorPlot = function(rsf.obj, 
                     eval.obj,
                     xlims,
                     ylims,
                     dir.results){
  
  # Plot 
  pdf(file.path(dir.results, paste0(rsf.obj$model, "_prediction_error_plot.pdf")))
  plot(eval.obj$pred.err, xlim =  xlims , ylim = ylims)
  dev.off()
}


# plotKMbyRiskGroup -------------------------------------------------------------------------


#
# Helper function of time conversion
# Input variable is days
timeConversion = function(x, to = "months") {
  if (to == "months") {
    return(x/30)
  } else {
    return(x/365)
  }
}

#
# Function will group individuals based on median predicted risk (mortality)
# into two groups (LOW/HIGH risk) and produce a Kaplan-Meier-plot
#
plotKMbyRiskGroup = function(rfs.obj,
                             end.point.event,
                             end.point.time,
                             time = "days",
                             plot.title){
  
  # Predicted risk (mortality)
  predicted.risk = rfs.obj$predicted
  
  # Find median risk 
  median.risk = median(predicted.risk)
  
  # Categorise into two groups 
  risk.group = ifelse(predicted.risk > median.risk, "High", "Low")
  
  # Survival data 
  survival.data = rfs.obj$yvar
  rownames(survival.data) = NULL
  
  # Create data frame 
  input.data = cbind(data.frame(risk.group = risk.group), survival.data)
  
  # Do time conversion if requested 
  if (time == "months") {
    input.data[end.point.time] = timeConversion(input.data[end.point.time], to = "months")
  } else if (time == "years") {
    input.data[end.point.time] = timeConversion(input.data[end.point.time], to = "years")
  } else {
    # Nothing
  }
  
  # Generate the formula for the model and fit 
  surv.expression = paste0("Surv(", end.point.time, ", " , end.point.event, ")")
  f <- as.formula(paste(surv.expression, "risk.group", sep = " ~ "))
  s.fit <- surv_fit(f, data = input.data)
  
  
  # Prepare KM
  s.fit.plot <- ggsurvplot(s.fit, 
                           data = input.data, legend = "bottom",
                           title = paste0(plot.title, " (n = ", nrow(input.data) ,")"),
                           legend.title = "Predicted risk", pval = TRUE, xlab = "Time (days)", 
                           font.family = "Helvetica", font.x = 15, font.y = 15, font.tickslab = 15, font.legend = 15,
                           conf.int = T,
                           break.time.by = 250,
                           surv.plot.height = 0.65, risk.table = TRUE, cumevents = F,
                           ggtheme = theme_classic(), 
                           fontsize = 5, pval.size = 7, tables.font.tickslab = 5, tables.y.text.col = T, tables.y.text = FALSE, 
                           tables.theme = theme_survminer())
  return(s.fit.plot)
}


# pltROC -------------------------------------------------------------------------

#
# Plot ROC curve 
#
#
pltROC = function(rsf.obj, 
                  validation.data,
                  end.point.event = "PFI.clin",
                  end.point.time = "PFI.time.clin",
                  fixed.times = NULL,
                  seed){
  
  # Set seed
  set.seed(seed)

  # Generate formula
  features = colnames(train_and_validation$validation)
  features_pred = features[!features %in% c("PFI.time.clin", "PFI.clin")]
  form = as.formula(paste(paste0("Surv(", end.point.time, "," , end.point.event, ")"), paste(features_pred, collapse = "+"), sep = "~"))
  
  # Calculate ROC scores 
  rR.score = riskRegression::Score(rsf.obj$fit.model, 
                                   formula = form, 
                                   data =  validation.data, 
                                   cens.model = "marginal",
                                   metrics = "auc",  
                                   times = fixed.times,
                                   summary = c("risks","IPA","riskQuantile","ibs"), 
                                   conf.int = TRUE, plots = "ROC", seed = seed)
  
  
  # Prepare plots 
  plot.ls = list()
  i = 1
  for (time in fixed.times) {
    
    # Extract data for plotting 
    roc.data = rR.score[["ROC"]]$plotframe %>% 
      dplyr::filter(times == time)
    
    # Get AUC and confidence interval 
    auc.data = rR.score$AUC$score %>% 
      dplyr::filter(times == time)
    
    auc = paste0(round(auc.data$AUC,3), " [", round(auc.data$lower,3)  ,"-", round(auc.data$upper,3) ,"]")
    
    # Prep plot
    gg = ggplot(roc.data, aes(x = 100 * FPR, y = 100 * TPR)) + geom_line(color = "red") + 
      geom_segment(aes(x=0,xend=100,y=0,yend=100), color = "skyblue", linetype = 2) +
      theme_classic() + xlab("1 - Specificity") + ylab("Sensitivity") + 
      annotate(geom = "text", 
               label = paste0("AUC = ", auc), 
               x = 20, y = 90)
    
    plot.ls[[i]] = gg
    i = i + 1
  }
  names(plot.ls) = paste0("t", fixed.times)
  
  
  return(plot.ls)
}



# prepareHeatmap ----------------------------------------------------------

# Helper function
recodeEvent = function(x){
  return(ifelse(x == 0, "Censored", "Progressed"))
}

# Convert numeric CN status to more descriptive 
convertCN = function(x){
  # Conversion of 
  conversion.table = list("-1" = "Deletion",
                          "0" = "Neutral" ,
                          "1" = "Amplification")
  return(unlist(conversion.table[x]))
}

# Helper function to cluster CN data
spearmanDist = function(x, y){
  # Conversion of 
  conversion.table = list("Deletion" = 1,
                          "Neutral" = 2,
                          "Amplification" = 3)
  
  x = unlist(conversion.table[x])
  y = unlist(conversion.table[y])
  return(sum((x - y) ** 2))
}

# Main function
prepareHeatmap = function(rsf.obj, dir.results = result.dir, prefix){
  
  # Feature data 
  feature.data = rsf.obj$xvar
  
  # Get CN data, convert to char and to matrix 
  feature.data.cn = dplyr::select(feature.data, ends_with(".cn"))
  cn.features.present = FALSE
  if (base::ncol(feature.data.cn) != 0) {
  	cn.mat = t(as.matrix(feature.data.cn))
  	cn.char = as.character(cn.mat )
  	cn.mat.discrete = matrix(cn.char, nrow = nrow(cn.mat))
  	cn.mat.discrete = apply(cn.mat.discrete, 1:2, convertCN)
  	colnames(cn.mat.discrete ) = NULL
  	row.names(cn.mat.discrete ) = str_remove(row.names(cn.mat), ".cn")

  	# Discrete colorscale for CN data
  	cn.col = structure(c("red","white","blue"),
                     names = c("Amplification", "Neutral", "Deletion"))

	cn.features.present = TRUE
  } 
  
  # Expression data, scale and convert to matrix
  feature.data.exp = dplyr::select(feature.data, ends_with(".exp"))
  exp.features.present = FALSE
  if (base::ncol(feature.data.cn) != 0) {
  	expression.mat = t(scale(as.matrix(feature.data.exp))) # Row scale 
  	colnames(expression.mat) = NULL
  	row.names(expression.mat) = str_remove(row.names(expression.mat), ".exp")
  
  	# Continuous color scale for expression
  	exp.col.fun = colorRamp2(seq(min(expression.mat), max(expression.mat), 
                               length = 3), 
                           c("blue", "white", "red"), space = "RGB")
  
  	exp.features.present = TRUE
  }


  print(exp.features.present)
  print(cn.features.present)
  # Extract survival data
  survival.data = rsf.obj$yvar
  
  # Recode the numberic data into factor
  survival.data = mutate(survival.data, Event = recodeEvent(PFI.clin))
  
  # Extract predicted risk 
  predicted.risk  = rsf.obj$predicted
  
  # Color scale for the prediction
  pred.col = colorRamp2(c(min(predicted.risk),max(predicted.risk)), 
                        c("#00bfc4","#f8766d"))
  
  # Calculate median risk
  median.risk = median(predicted.risk)
  # Categorise the predicted risk into High and Low based on the median
  risk.group = ifelse(predicted.risk >= median.risk, "High", "Low")
  # Colorscale for the risk groups
  risk.col = structure(c("#00bfc4","#f8766d"), names = c("Low", "High"))
  
  # Colorscale for the time variable
  time.col = colorRamp2(seq(0, max(survival.data$PFI.time.clin), length = 2), 
                        c("white", "green"), space = "RGB")
  
  # Colorscale for the status 
  status.col = structure(c("grey","black"), names = c("Censored","Progressed"))
  
  # Prepare the heatmap components  
  
  # Top annotation
  ha = HeatmapAnnotation(`Predicted risk` = anno_barplot(predicted.risk, border = F, axis = F, bar_width = 1,
                                                         gp = gpar(fill = 4, lwd = 0)),
                         `Time (days)` = survival.data$PFI.time.clin,
                         Status = survival.data$Event,
                         col = list(Status = status.col,
                                    `Time (days)` = time.col))
  
  # Bottom annotation
  ha2 = HeatmapAnnotation(`Risk group` = risk.group,
                          col = list(`Risk group` = risk.col))
  
  if ((cn.features.present == T) & (exp.features.present == T)) {
	
       expression.heat = Heatmap(expression.mat,
                            name = "Expression",
                            row_title = "Expression",
                            column_order = order(predicted.risk),
                            top_annotation = ha,
                            col = exp.col.fun,
                            border = T,
                            row_names_gp = gpar(fontsize = 8),
                            height = nrow(expression.mat)*unit(4, "mm"))

	cn.heat = Heatmap(cn.mat.discrete,
                    name = "CN status",
                    row_title = "CN status",
                    column_order = order(predicted.risk),
                    clustering_distance_rows = spearmanDist,
                    col = cn.col,
                    bottom_annotation = ha2,
                    border = T,
                    row_names_gp = gpar(fontsize = 8),
                    height = nrow(cn.mat.discrete)*unit(4, "mm"))


  	# Combine the heatmaps 
  	comb.heat = expression.heat %v% cn.heat


  } else if ((cn.features.present == F) & (exp.features.present == T)) {
  	
	expression.heat = Heatmap(expression.mat,
                            name = "Expression",
                            row_title = "Expression",
                            column_order = order(predicted.risk),
                            top_annotation = ha,
                            bottom_annotation = ha2,
                            col = exp.col.fun,
                            border = T,
                            row_names_gp = gpar(fontsize = 8),
                            height = nrow(expression.mat)*unit(4, "mm"))


  	# Combine the heatmaps 
  	comb.heat = expression.heat

  } else {
	
	cn.heat = Heatmap(cn.mat.discrete,
                    name = "CN status",
                    row_title = "CN status",
                    column_order = order(predicted.risk),
                    clustering_distance_rows = spearmanDist,
                    col = cn.col,
                    top_annotation = ha,
                    bottom_annotation = ha2,
                    border = T,
                    row_names_gp = gpar(fontsize = 8),
                    height = nrow(cn.mat.discrete)*unit(4, "mm"))


  	# Combine the heatmaps 
  	comb.heat = cn.heat

  }
  
  pdf(file.path(dir.results, paste0(prefix, "_heatmap.pdf")))
  draw(comb.heat)
  dev.off()
}

# plotIFE --------------------------------------------------------------------
#
# Function for visualisation of results from the iterative feature elimination
#

# Prepare visualisation
plotIFE = function(ife.results, dir.results) {
   
    # Prepare plot
    gg = ggplot(ife.results, 
	    aes(x = num_variable, y = error)) + 
            geom_line() + theme_classic() + 
	    ylim(0, 1)

    ggsave(file.path(dir.results, "IFE_results.pdf"))
}




