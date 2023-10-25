# Scripts for plotting a complex heatmap of the glmnet results for TCGA-KIRC
# Code adapted from functions by Tommi Rantapero 

# Helper function to recode event variable 
recodeEvent = function(x){
  return(ifelse(x == 0, "Censored", "Progressed"))
}

# Convert numeric CN status to more descriptive 
convertCN = function(x){
  # Conversion of 
  conversionTable = list("-1" = "Deletion",
                          "0" = "Neutral" ,
                          "1" = "Amplification")
  return(unlist(conversionTable[x]))
}

# Helper function to cluster CN data
spearmanDist = function(x, y){
  # Conversion of 
  conversionTable = list("Deletion" = 1,
                          "Neutral" = 2,
                          "Amplification" = 3)
  
  x = unlist(conversionTable[x])
  y = unlist(conversionTable[y])
  return(sum((x - y) ** 2))
}

#
# Main function for preparing heatmap
#
# UPDATED VERSION 27.9.2022
#
prepareHeatmap = function(feature_data, 
                          surv_data, 
                          predicted.risk, 
                          dir, 
                          name, 
                          row.height = 4){

  # Get CN data, convert to char and to matrix 
  featureDataCN = dplyr::select(feature_data, ends_with(".cn"))
  cn.features.present = FALSE
  if (base::ncol(featureDataCN) != 0) {
    cnMat = t(as.matrix(featureDataCN))
    cnChar = as.character(cnMat)
    cnMatDiscrete = matrix(cnChar, nrow = nrow(cnMat))
    cnMatDiscrete = apply(cnMatDiscrete, 1:2, convertCN)
    colnames(cnMatDiscrete ) = NULL
    row.names(cnMatDiscrete ) = str_remove(row.names(cnMat), ".cn")
    
    # Discrete colorscale for CN data
    cnCol = structure(c("red","white","blue"),
                       names = c("Amplification", "Neutral", "Deletion"))
    
    cn.features.present = TRUE
  } 
  
  # Expression data, scale and convert to matrix
  featureDataExp = dplyr::select(feature_data, ends_with(".exp"))
  exp.features.present = FALSE
  if (base::ncol(featureDataExp) != 0) {
    expressionMat = t(scale(as.matrix(featureDataExp))) # Row scale 
    colnames(expressionMat) = NULL
    row.names(expressionMat) = str_remove(row.names(expressionMat), ".exp")
    
    # Continuous color scale for expression
    expColFun = colorRamp2(seq(min(expressionMat), max(expressionMat), 
                                 length = 3), 
                             c("blue", "white", "red"), space = "RGB")
    
    exp.features.present = TRUE
  }
    
  if ((cn.features.present == T) | (exp.features.present == T)) {
      
        # Extract survival data
        survivalData = surv_data
  
        # Recode the numberic data into factor
        survivalData <- mutate(survivalData, Event = recodeEvent(status))

        # Extract predicted risk 
        predRisk  = predicted.risk

        # Based on the predicted risk divide the patients into groups 
        group = ifelse(predRisk  < median(predRisk), "Low", "High")
  
        # Color scale for the prediction
        predCol = colorRamp2(c(min(predRisk),max(predRisk)), 
                        c("#00bfc4","#f8766d"))

        # Colorscale for the risk groups
        riskCol = structure(c("#00bfc4","#f8766d"), names = c("Low", "High"))


        # Colorscale for the time variable
        timeCol = colorRamp2(seq(0, max(survivalData$time), length = 2), 
                        c("white", "green"), space = "RGB")
  
        # Colorscale for the status 
        statusCol = structure(c("grey","black"), names = c("Censored","Progressed"))

        # Initialise clinical variables and colorschemes 
        age.at.diag.years = NULL
        gender = NULL
        tumor.stage = NULL
        gleason.score = NULL 
    
        ageCol = NULL
        genderCol = NULL
        tumorStageCol = NULL 
        gleasonCol = NULL

        # Add clinical variables if present 
        if (is.null(feature_data$Age) == F) {
            age.at.diag.years = (feature_data$Age / 360) 
        }

        gender = feature_data$Gender
        tumor.stage = feature_data$Tumor.stage
        gleason.score = feature_data$Gleason.group
 
        # Age colorscheme
        if (is.null(feature_data$Age) == F) {
            ageCol = colorRamp2(c(min(age.at.diag.years),max(age.at.diag.years)), 
                        c("white","midnightblue"))
        }

        # Gender colorscheme
        if (is.null(feature_data$Gender) == F){
            genderCol = structure(c("orange", "brown"), names = c("female","male"))
        }

        # Tumor stage colorscheme
        if (is.null(feature_data$Tumor.stage) == F){
            tumorStageCol = structure(c("grey","steelblue1", "steelblue2",  "steelblue3", "steelblue"),
                                   names = c("NA","Stage 1", "Stage 2", "Stage 3", "Stage 4"))}

        if (is.null(feature_data$Gleason.group) == F) {
            gleasonCol = structure(c("seagreen3","tomato2"), names = c("Gleason group 1","Gleason group 2"))
        }

        # Prepare the heatmap components  

        # Collect varibles into a final data frame 
        heatmap.data = data.frame(`Time (days)` = survivalData$time,
                              Status = survivalData$Event, check.names = F)

        if (is.null(age.at.diag.years) == F) {
            heatmap.data$`Age at diagnosis` = age.at.diag.years
        }
        if (is.null(gender) == F) {
            heatmap.data$`Gender` = gender
        }
        if (is.null(tumor.stage) == F) {
            heatmap.data$`Tumor stage` = tumor.stage
        } 

        if (is.null(gleason.score) == F) {
            heatmap.data$`Gleason score` = gleason.score
        } 

        # Color data for the heatmap
        color.data = list(Status = statusCol,
                      `Time (days)` = timeCol)

        if (is.null(age.at.diag.years) == F) {
            color.data$`Age at diagnosis` = ageCol
        }
        if (is.null(gender) == F) {
            color.data$Gender = genderCol
        }
        if (is.null(tumor.stage) == F) {
            color.data$`Tumor stage` = tumorStageCol
        } 
        if (is.null(gleason.score) == F) {
            color.data$`Gleason score` = gleasonCol
        } 

        # Generate top annotation 
        ha.top.1 = HeatmapAnnotation(`Predicted risk` = anno_barplot(predRisk, 
                                                                 border = F, 
                                                                 axis = F, 
                                                                 bar_width = 1,
                                                  gp = gpar(fill = 4, lwd = 0)))

        ha.top.2 = HeatmapAnnotation(df = heatmap.data, 
                           col = color.data)

        # Concatenate 
        ha.top = c(ha.top.1, ha.top.2) 

        # Bottom annotation
        ha.bottom = HeatmapAnnotation(`Risk group` = group,
                          col = list(`Risk group` = riskCol))
    
        if ((cn.features.present == T) & (exp.features.present == T)) {
    
            expressionHeat = Heatmap(expressionMat,
                              name = "Expression",
                              row_title = "Expression",
                              column_order = order(predRisk),
                              top_annotation = ha.top,
                              col = expColFun,
                              border = T,
                              row_names_gp = gpar(fontsize = 8),
                              height = nrow(expressionMat)*unit(row.height, "mm"))
    
            cnHeat = Heatmap(cnMatDiscrete,
                      name = "CN status",
                      row_title = "CN status",
                      column_order = order(predRisk),
                      clustering_distance_rows = spearmanDist,
                      col = cnCol,
                      bottom_annotation = ha.bottom,
                      border = T,
                      row_names_gp = gpar(fontsize = 8),
                      height = nrow(cnMatDiscrete)*unit(row.height, "mm"))
    
    
            # Combine the heatmaps 
            combHeat = expressionHeat %v% cnHeat
            
            pdf(file.path(dir, paste0(name, "_heatmap.pdf")))
            draw(combHeat)
            dev.off()
            return(combHeat)
    
    
        } else if ((cn.features.present == F) & (exp.features.present == T)) {
    
            expressionHeat = Heatmap(expressionMat,
                              name = "Expression",
                              row_title = "Expression",
                              column_order = order(predRisk),
                              top_annotation = ha.top,
                              bottom_annotation = ha.bottom,
                              col = expColFun,
                              border = T,
                              row_names_gp = gpar(fontsize = 8),
                              height = nrow(expressionMat)*unit(row.height, "mm"))
    
    
            # Combine the heatmaps 
            combHeat = expressionHeat
        
            pdf(file.path(dir, paste0(name, "_heatmap.pdf")))
            draw(combHeat)
            dev.off()
            return(combHeat)

        } else {
            
            cnHeat = Heatmap(cnMatDiscrete,
                      name = "CN status",
                      row_title = "CN status",
                      column_order = order(predRisk),
                      clustering_distance_rows = spearmanDist,
                      col = cnCol,
                      bottom_annotation = ha.bottom,
                      border = T,
                      row_names_gp = gpar(fontsize = 8),
                      height = nrow(cnMatDiscrete)*unit(row.height, "mm"))

            combHeat = cnHeat
        
            pdf(file.path(dir, paste0(name, "_heatmap.pdf")))
            draw(combHeat)
            dev.off()
            return(combHeat)

        }
      
    } else {
        return(NULL)
    }
}





