#
# TCGA data preprocessing functions
#
#

library("DESeq2")

#
# Subset TCGA object. This is supposed to be a general 
# utility function to remove genes or samples from 
# the dataset. 
#
#

subsetTCGAdataObj = function(data.obj, 
                             by.gene.id = NULL,  
                             by.gene.symbol = NULL,  
                             by.sample = NULL) {
  
  # If by.gene.id is not null subset 
  if (is.null(by.gene.id) == F) {
    
    # Which rows to keep
    rows.keep = data.obj$Annotation$Gene.ID %in% by.gene.id
    
    # This will only affect the Data and Annotation elements 
    data.obj$Annotation = data.obj$Annotation[rows.keep, ]
    data.obj$Data = data.obj$Data[rows.keep,]
    
  } 
  
  # If by.gene.id is not null subset 
  if (is.null(by.gene.symbol) == F) {
    
    # Which rows to keep
    rows.keep = data.obj$Annotation$Gene.Symbol %in% by.gene.symbol
    
    # This will only affect the Data and Annotation elements 
    data.obj$Annotation = data.obj$Annotation[rows.keep, ]
    data.obj$Data = data.obj$Data[rows.keep,]
    
  } 
  
  
  # If by.sample.id is not null subset 
  if (is.null(by.sample) == F) {
    
    # Which columns to keep
    cols.keep = data.obj$Data %in% by.sample
    
    # This will affect the Data. Technically it 
    # can affect clinical data but since patient can have 
    # many samples we will deal with this later
    data.obj$Data = data.obj$Data[,cols.keep]
    
  }
  return(data.obj)
}

#
# Drop out normal samples  
# The sample codes for normal samples are :
# 10, 11, 12, 13, 14
dropNormalSamples = function(data.obj) {
  
  # Helper function for getting the sample type code
  getSampleType = function(sample.id){
    # Extract sample code
    sample.code = unlist(stringr::str_split(string = sample.id, 
                                            pattern = "-"))[4]
    # Get numeric code 
    sample.code.num = stringr::str_extract(sample.code, "\\d+")
    
    return(sample.code.num)
    
  }
  
  # Get the sample type codes 
  sample.codes = map_chr(colnames(data.obj$Data), getSampleType)

  
  # Helper function to test if sample is derived from normal tissue
  isNormal = function(x) {
    if (x %in% c("10","11","12","13","14")) {
      return(T)
    } else {
      return(F)
    }
  }
  
  # Samples to keep
  keep.sample = !(unlist(map(sample.codes, isNormal)))
  
  
  
  # Subset the Data 
  data.obj$Data = data.obj$Data[,keep.sample]
  
  # If Raw slot present we filter also this data
  if ("Raw" %in% names(data.obj)) {
    data.obj$Raw = data.obj$Raw[,keep.sample]
  }
  
  return(data.obj)
}

#
# Function for selecting only the primary tumor samples
#
#

# Function for selecting only primiary tumor samples
selectPrimaryT = function(data.obj){
  
  # Helper function for getting the sample type code
  getSampleType = function(sample.id){
    # Extract sample code
    sample.code = unlist(stringr::str_split(string = sample.id, 
                                            pattern = "-"))[4]
    # Get numeric code 
    sample.code.num = stringr::str_extract(sample.code, "\\d+")
    
    return(sample.code.num)
    
  }
  
  # Get the sample type codes 
  sample.codes = map_chr(colnames(data.obj$Data), getSampleType)

  # Helper function to test if sample is derived from normal tissue
  isPrimary = function(x) {
    if (x %in% c("01")) {
      return(T)
    } else {
      return(F)
    }
  }
  
  # Samples to keep
  keep.sample = unlist(map(sample.codes, isPrimary))
  
  # Subset the Data 
  data.obj$Data = data.obj$Data[,keep.sample]
  
  # If Raw slot present we filter also this data
  if ("Raw" %in% names(data.obj)) {
    data.obj$Raw = data.obj$Raw[,keep.sample]
  }
  
  return(data.obj)
}

#
# Special function for Acute myoloid leukemia
# THere are no primary solid tumor so we 
# select samples with code 03 = Primary Blood Derived Cancer
selectPrimaryB = function(data.obj){

  # Helper function for getting the sample type code
  getSampleType = function(sample.id){
    # Extract sample code
    sample.code = unlist(stringr::str_split(string = sample.id, 
                                            pattern = "-"))[4]
    # Get numeric code 
    sample.code.num = stringr::str_extract(sample.code, "\\d+")
    
    return(sample.code.num)
    
  }
  
  # Get the sample type codes 
  sample.codes = map_chr(colnames(data.obj$Data), getSampleType)

  # Helper function to test if sample is derived from normal tissue
  isPrimary = function(x) {
    if (x %in% c("03")) {
      return(T)
    } else {
      return(F)
    }
  }
  
  # Samples to keep
  keep.sample = unlist(map(sample.codes, isPrimary))
  
  # Subset the Data 
  data.obj$Data = data.obj$Data[,keep.sample]
  
  # If Raw slot present we filter also this data
  if ("Raw" %in% names(data.obj)) {
    data.obj$Raw = data.obj$Raw[,keep.sample]
  }
  
  return(data.obj)

}



#
# Drop out duplicated samples
# The sample ids are converted into participant ids 
# In the case of multiple samples we prioritise 01A sample  
# If there are multiple such samples we select the first one we 
# encounter
#

dropDuplicateSamples = function(data.obj) {

  # Get the sample ids 
  sample.ids = colnames(data.obj$Data)
  
  # Find the duplicated entries 
  participant.ids = unlist(map(sample.ids, sampleToParticipantID))
  
  # Counts 
  duplicated.ids = names(table(participant.ids)[table(participant.ids) > 1])
  non.duplicated.ids = names(table(participant.ids)[table(participant.ids) == 1])

  # Return the indices of the duplicated samples 
  dup.idx = which(participant.ids %in% duplicated.ids)
  non.dup.idx = which(participant.ids %in% non.duplicated.ids)

  # Sample ids of the duplicated samples 
  duplicated.samples = sample.ids[dup.idx]
  non.duplicated.samples = sample.ids[non.dup.idx]

  # Find 01A samples 
  samples.01A = duplicated.samples[str_detect(duplicated.samples, "01A")]

  # Still from this set remove the duplicates 
  samples.01A.participant.id = unlist(map(samples.01A, sampleToParticipantID))
   
  # Pick the first instance 
  pickFirstInstance = function(part.id, sample.ids) {
    return(sample.ids[which(str_detect(sample.ids, part.id))[1]])
  }

  samples.keep = map_chr(unique(samples.01A.participant.id), pickFirstInstance, sample.ids = samples.01A)
  samples.keep = c(samples.keep, non.duplicated.samples)

  # Remove samples from clinical and CN data 
  data.obj$Clinical.data = data.obj$Clinical.data[rownames(data.obj$Clinical.data) %in% samples.keep,]
  data.obj$Data = data.obj$Data[,colnames(data.obj$Data) %in% samples.keep]
  # If Raw slot present filter also 
  if ("Raw" %in% names(data.obj)) {
    data.obj$Raw = data.obj$Raw[,colnames(data.obj$Raw) %in% samples.keep]
  }
  return(data.obj)
}


#
# UPDATE 22.1.2024 
#
#
# Based on the reviewer comments we will do VST separately 
# for training and test sets according to :
# https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation
#

#
# Perform variance stabilising transformation for 
# gene expression data for the Training set 
#

performVSTtraining = function(counts){

  # Prepare minimal colData with only one condition
  coldata = data.frame("Condition" = rep("Condition", ncol(counts)))
  rownames(coldata) = colnames(counts)

  # Generate 
  dds = DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~1)

  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)

  # Variance stabilisation function 
  vds  = vst(dds, blind = T)

  # Get the dispersion function to be 
  # used with also the test set 
  disp.function = dispersionFunction(dds)
  
  # Get the transformed counts
  vst.counts = assay(vds)
 
  # Prepare and object with the vst counts 
  # for the training set and the estimated dispersion function 
  vst.object = list(vst.counts = vst.counts, 
                    disp.function = disp.function) 
  
  return(vst.object)
   
}

performVSTtest = function(counts, disp.function.train){

  # Prepare minimal colData with only one condition
  coldata = data.frame("Condition" = rep("Condition", ncol(counts)))
  rownames(coldata) = colnames(counts)

  # Generate 
  dds = DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~1)

  dds = estimateSizeFactors(dds)
  #dds = estimateDispersions(dds)
  dispersionFunction(dds) = disp.function.train 

  # Variance stabilisation function 
  vds  = vst(dds, blind = T)

  # Get the transformed counts
  vst.counts = assay(vds)
 
  return(vst.counts)
   
}

#
# Old function for performing the VST for full dataset
#
performVST = function(counts){
 
   # Prepare minimal colData with only one condition
   coldata = data.frame("Condition" = rep("Condition", ncol(counts)))
   rownames(coldata) = colnames(counts)

   # Generate
   dds = DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                colData = coldata,
                                design = ~1)

   # Variance stabilisation function
   vds  = vst(dds, blind = T)

   # Get the transformed counts
   vst.counts = assay(vds)

   return(vst.counts)
}

#
# Function for preparing data matrix 
# from CN data. The idea is that the gene 
# symbols are added as row names. Possible duplicates
# or entries without symbol are treated the following way
# No gene symbol => ensembl id
# Duplicate gene symbol => <gene.symbol>.rowindex
#
#
prepareDataMatrix = function(tcga.data) {
  
  # Get the Data
  data.df = as.data.frame(tcga.data$Data)
  
  # Prepare row names. We will use the gene symbol if available. Otherwise we 
  # assign anon + row name. This way we can obtain some information 
  
  # Helper function
  prepRowNames = function(symbol, ens.id){
    if ((is.na(symbol) != T) & (symbol != "NA")){
      return(symbol)
    } else {
      return(ens.id)
    }
  }
  
  # Get gene names 
  rownames = unlist(purrr::map2(tcga.data$Annotation$Gene.Symbol,tcga.data$Annotation$Gene.ID, prepRowNames))
  
  # If there is duplicate gene ids add row index 
  rownames.count = table(rownames)
  
  # Helper function to resolve duplicates
  resolveDuplicates = function(rowname, idx ,counts){
    if (counts[[rowname]] > 1) {
      return(paste0(rowname, ".", idx))
    } else {
      return(rowname)
    }
  } 
  rownames.uniq = unlist(purrr::imap(rownames, resolveDuplicates, counts =  rownames.count))
  rownames(data.df) =  rownames.uniq
  
  return(data.df)
}

#
# The same exact function as above but for raw data 
#
prepareDataMatrixRaw = function(tcga.data) {
  
  # Get the Data
  data.df = as.data.frame(tcga.data$Raw)
  
  # Prepare row names. We will use the gene symbol if available. Otherwise we 
  # assign anon + row name. This way we can obtain some information 
  
  # Helper function
  prepRowNames = function(symbol, ens.id){
    if ((is.na(symbol) != T) & (symbol != "NA")){
      return(symbol)
    } else {
      return(ens.id)
    }
  }
  
  # Get gene names 
  rownames = unlist(purrr::map2(tcga.data$Annotation$Gene.Symbol, tcga.data$Annotation$Gene.ID, prepRowNames))
  
  # If there is duplicate gene ids add row index 
  rownames.count = table(rownames)
  
  # Helper function to resolve duplicates
  resolveDuplicates = function(rowname, idx ,counts){
    if (counts[[rowname]] > 1) {
      return(paste0(rowname, ".", idx))
    } else {
      return(rowname)
    }
  } 
  rownames.uniq = unlist(purrr::imap(rownames, resolveDuplicates, counts =  rownames.count))
  rownames(data.df) =  rownames.uniq
  
  return(data.df)
}



#
# Helper findividual. Extact the sample id from the aliquot id 
#
aliquotToSampleID = function(data.mat){
  
  # Sample ids 
  sample.ids = colnames(data.mat)
  
  # Helper function
  getID = function(sample.id){
    id.split = unlist(stringr::str_split(string = sample.id, pattern = "-"))
    return(paste(id.split[1:4], collapse = "-"))
  }

  partID = unlist(map(sample.ids, getID))
  return(partID)
  
}


#
# Helper function. Extract the participant id from the sample id. This 
# basically defines the individual
#
sampleToParticipantID = function(sample.id){
  id.split = unlist(stringr::str_split(string = sample.id, pattern = "-"))
  part.id = paste(id.split[1:3], collapse = "-")
  return(part.id)
}  
    

#
# Function for merging experimental data and clinical data 
# 

mergeTCGAdata = function(clin.data, data.matrices, data.suffixes){
  
  # Clinical data 
  colnames(clin.data) = paste0(colnames(clin.data), ".clin") 
  
  # Convert bcr_patient_barcode.clin to Participant.ID
  clin.data = dplyr::rename(clin.data, "Participant.ID" = "bcr_patient_barcode.clin") 
  
  # Add suffixes to data matrices
  for (i in 1:length(data.matrices)) {
    rownames(data.matrices[[i]]) = paste0(rownames(data.matrices[[i]]), ".", data.suffixes[i]) 
  }
  
  # Convert the data matrix column names to sample names 
  for (i in 1:length(data.matrices)) {
    colnames(data.matrices[[i]]) = aliquotToSampleID(data.matrices[[i]])
  }
  
  # Merge the experimental data 
  merged.data = purrr::reduce(data.matrices, dplyr::bind_rows)
  
  # Transpse 
  merged.data = t(merged.data) %>% 
                     as.data.frame() %>% 
                     tibble::rownames_to_column(var = "Sample.ID")
  
  # Add the participant ID 
  merged.data["Participant.ID"] = map_chr(merged.data$Sample.ID, sampleToParticipantID)
  
  # Merge clinical data and 
  merged.data.final = dplyr::left_join(clin.data, merged.data, by = "Participant.ID")
  
    
  # Split again based on the data type and generate the final object 
  merged.obj = list("CLIN" =  merged.data.final %>% dplyr::select(Participant.ID, contains(".clin")), 
                    "EXP" = merged.data.final %>% dplyr::select(Participant.ID, contains(".exp")),
                    "CN" = merged.data.final %>% dplyr::select(Participant.ID, contains(".cn")))   
  
  return(merged.obj)
}

# --- FUNCTION FOR SELECTION OF PATIENTS AND VARIABLES 

selectVariables = function(clinical.endpoint, 
                           gene.list, 
                           data.suffixes){

    # Construct the clinical features   
    clinical.features =  c(paste0(clinical.endpoint, ".clin"),
                           paste0(clinical.endpoint, ".time.clin"))

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
splitCases = function(data, split, only.complete = T, variables, seed, impute.valid = F, impute.var = NULL) {

  
  # Number of cases in total
  ncases.initial = nrow(data)
  
  if (only.complete == T) {
     # Take only complete cases if data contains missing values 
     print("Taking only complete cases")
        
     variables = variables[variables %in% colnames(data)]
     complete <- complete.cases(data[,variables]) # Indexes of cases with no missing values among input variables
     data <- data[complete, c("submitter_id.clin", variables)]
     rownames(data) = NULL
     data = tibble::column_to_rownames(data, var = "submitter_id.clin")
     print(paste0("Including ", sum(complete)  ," cases out of ", ncases.initial, " cases"))
  } else {
     # Allow missing values  	  
     print("Taking all cases")
     data <- data[,c("submitter_id.clin", variables)]
     rownames(data) = NULL
     data = tibble::column_to_rownames(data, var = "submitter_id.clin")
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


# --- FUNCTIONS FOR PROCESSING VARIABLES --- #

#
# Reformat variable Tumor Stage : ajcc_pathologic_stage 
#
 
truncateTumorStage = function(stage){
    if (is.na(stage) == T){
        return("NA")
    } else {
        return(stringr::str_extract(stage, "T\\d"))
    }
}

getTumorStage = function(x){

    # Truncate the T 
    trunc.stage = map_chr(x$ajcc_pathologic_t.clin, truncateTumorStage)

    # Convert into factor 
    tumor.stage = factor(trunc.stage, levels = c("T0","T1","T2","T3","T4"))

    # Return
    return(tumor.stage)
}

# Set up variables of interest with appropriate factor levels 
reformatTumorStage = function(x){

    # Patterns to match
    pattern = c("[IV]+")

    # Extract pattern 
    stage = str_extract(x, pattern)

    if ((is.na(stage) == T) | (stage == "Stage X") | (stage == "Stage Tis") | (stage == "Unknown") | (stage == "Not Reported")) { 
        return(NA)
    }

    if (stage == "I") {
        return("Stage 1")
    } else if (stage == "II") {
        return("Stage 2")
    } else if (stage == "III") {
        return("Stage 3")
    } else {
        return("Stage 4")
    }
}

#
# Find the Gleason group 
# 

determineGleasonGroup = function(primary.patt, secondary.patt){
    primary.score = as.numeric(stringr::str_extract(primary.patt, "\\d"))
    secondary.score = as.numeric(stringr::str_extract(secondary.patt, "\\d"))
    gleason.score = primary.score + secondary.score
    if (gleason.score < 7){
       return("Gleason_group_1") 
    } else if (gleason.score > 7) {
       return("Gleason_group_2") 
    } else {
        if (primary.score == 4) {
            return("Gleason_group_2")
        } else {
            return("Gleason_group_1")
        }
    }
}

#
# Convert Copy-number status variables to factors 
#

# Convert cn-features to factors 
convertCnToFactors = function(x){
    # Change neutral CN as the reference level
    x = dplyr::mutate(x, across(contains(".cn"), 
                               function(x){factor(x, levels = c("0", "-1", "1"))})) 
   return(x)
}


# --- FUNCTIONS FOR CALCULATING STATS --- #


#
# Helper functions
#


# Counting the fraction of NAs 
countNaFrac = function(x){
    return(sum(is.na(x))/length(x))
}

# Counting the fraction of amplifications 
countAmpFrac = function(x){
   return(sum(x[is.na(x) == F]== 1)/length(x[is.na(x) == F]))
}

# Counting the fraction of deletions
countDelFrac = function(x){
   return(sum(x[is.na(x) == F] == -1)/length(x[is.na(x) == F]))
}


#
# Main function
#

prepSummaryCN = function(x, variables, type){

    # Select variables of interest 
    variables.cn = variables[str_ends(variables, ".cn")]
        
    # Check which feature are present 
    variables.cn = variables.cn[variables.cn %in% colnames(x)]

    # For expression data count the missing values, mean, and variance
    summary.table = x %>% 
                        dplyr::select(all_of(variables.cn)) %>% 
                        pivot_longer(everything()) %>%
                        group_by(name) %>% 
                        summarise(`Fraction with amplification` = countAmpFrac(value),
                                  `Fraction with deletion` = countDelFrac(value),
                                  `Fraction with NA` = countNaFrac(value))

    # Calculate the max 
    summary.table$`Maximum fraction of aberrations` = summary.table %>% 
                               dplyr::select(`Fraction with amplification`, `Fraction with deletion`) %>%
                               apply(1, max)  

    return(summary.table)
}

#
# Helper functions
#

countZeroFrac = function(x){
    res = sum(x[is.na(x) == F] == 0)/ length(x[is.na(x) == F])
    return(res)
}

#
# Main function
#

prepSummaryExp = function(x, raw.data, variables, type){

    # Select variables of interest 
    variables.exp = variables[str_ends(variables, ".exp")]

    # Check which feature are present 
    variables.exp = variables.exp[variables.exp %in% colnames(x)]

    # Modify the raw.data matrix 
    colnames(raw.data) = map_chr(colnames(raw.data), sampleToParticipantID)

    # Append the .exp extension 
    rownames(raw.data) = paste0(rownames(raw.data), ".exp")

    # Select only patients included in training set 
    raw.data = raw.data %>%
                        dplyr::select(one_of(tcga.dataset.merged[[end.point]]$train$Participant.ID))
    
    # Only include the selected variables 
    raw.data = raw.data[variables.exp,]

    # Calculate the statistics for the raw data 
    summary.table = raw.data %>%
        t() %>%
        as.data.frame() %>%
        pivot_longer(everything()) %>%
        group_by(name) %>% 
        summarise(`Median expression` = median(value, na.rm = T),
                  `Expression Variance` = var(value, na.rm = T),
                  `Fraction of zero expression`= countZeroFrac(value))

    return(summary.table)

}


#
# Functions for filtering features based on summary statistics
#

#
# Helper function
#
isInVect = function(x, type, feature.ls){
    if (str_ends(x, type) == T){
        if (x %in% feature.ls){
            return(T)
        } else {
            return(F)
        }
    } else {
        return(T)
    }
}

#
# MAIN function
#

# Keep only specified features by a given type
filterFeatures = function(x, keep, type){
    keep.feature = map_lgl(x, isInVect, type = type, feature.ls = keep)
    return(x[keep.feature])
}



# --- FUNCTIONS FOR PRODUCING INPUT FOR MODELLING --- #

#
# Helper function 
#

findReferenceLevel = function(x.train, x.validation){
    
    # Iterate over tumor stages and find the ones to be 
    # removed 
    rem.train = c()
    rem.valid = c()

    # Levels 
    levs = c("Stage 1","Stage 2","Stage 3","Stage 4")

    # Iterate over the tumor stages (Training data)
    for (stage in levs){
        col = paste0("Tumor.stage",stage)
        if (col %in% colnames(x.train)) {
            if (sum(x.train[,col]) == 0){
                rem.train = c(rem.train, col)
            } else {
                rem.train = c(rem.train, col)
                break
            }
        } else {
            rem.train = c(rem.train, col)
        }
    }
    
    # Iterate over the tumor stages (Validation data)
    for (stage in levs){
        col = paste0("Tumor.stage",stage)
        if (col %in% colnames(x.validation)) {
            if (sum(x.validation[,col]) == 0){
                rem.valid = c(rem.valid, col)
            } else {
                rem.valid = c(rem.valid, col)
                break
            }
        } else {
            rem.valid = c(rem.valid, col)
        }
    }

    # Find the intersection of the two 
    levs.rm = intersect(rem.train, rem.valid)
    return(levs.rm)
}


#
# Function for finding the reference level 
#

removeReferenceTumorStage = function(x, levs.rm){

    x = data.frame(x, check.names = F)
    
    # Levels 
    levs = paste0("Tumor.stage", c("Stage 1","Stage 2","Stage 3","Stage 4"))

    # We might possibly add zero columns for stages that are higher than 
    # reference if missing 
    for (lev in levs){
        if (!(lev %in% levs.rm)){
            if (!(lev %in% colnames(x))){
                print("jee")
                x[lev] = 0
            }
        }
    }
    
    # Finally remove the reference level and everyting lower
    x = x[, !(colnames(x) %in% levs.rm)]
    return(as.matrix(x))
}

# 
# Function removes the reference tumor stage. Notice that 
# there reference here is considered to be the lowest stage 
# which we have observed in the data. If there are any stages 
# lower than the reference stage those are removed too. This will
# leave only column for the higher tumor stages which are compared 
# against the reference
# 


#
# Function for generating model matrices
#

generateModelMatrices = function(train.data, 
                                 validation.data,
                                 clinical.endpoint){
    
    # Outcome variables
    time.var = paste0(clinical.endpoint, ".time.clin")
    status.var = paste0(clinical.endpoint, ".clin")
    
    # Remove possible zero time events which will cause errors.
    train.data <- dplyr::filter(train.data, !!rlang::sym(time.var) != 0)
    validation.data <- dplyr::filter(validation.data, !!rlang::sym(time.var) != 0)
    
    ### Outcome variables 
    
    # Extract events and time to event
    y.train = train.data %>% dplyr::select(all_of(c(status.var, time.var)))
    y.valid = validation.data %>% dplyr::select(all_of(c(status.var, time.var)))
    
    # Assign the column names to follow standard convetions
    colnames(y.train) = c("status", "time")
    colnames(y.valid) = c("status", "time")

    # Generate a survival object for fitting the model
    y.train <- Surv(y.train$time, y.train$status) 
    ###
    
    ### Predictors 
    
    # We keep the outcome variables in the model until the last step to
    # avoid collapsing the data into a vector given that only one predictor 
    # is present
    
    # Convert cn-features to factors 
    x.train = convertCnToFactors(train.data)
    x.valid = convertCnToFactors(validation.data)

    # Create a model matrix
    # The CN features will be converted to dummy variables as the glmnet does 
    # not deal with factors. NOTE: THERE MIGHT BE NEED TO CHANGE THIS ENCODING
    x.train.mat = model.matrix(~ .-1, x.train)
    x.valid.mat = model.matrix(~ .-1, x.valid)

    # Remove the Genderfemale column as it's the reference level for Gender
    x.train.mat = x.train.mat[, colnames(x.train.mat) != "Gender.clinfemale"]
    x.valid.mat = x.valid.mat[, colnames(x.valid.mat) != "Gender.clinfemale"]

    # Remove Gleason group 1 as it´s the reference level for Gleason group 
    x.train.mat = x.train.mat[, colnames(x.train.mat) != "Gleason.group.clinGleason_group_1"]
    x.valid.mat = x.valid.mat[, colnames(x.valid.mat) != "Gleason.group.clinGleason_group_1"]

    # Find the reference tumor stage 
    levs.rm = findReferenceLevel(x.train.mat, x.valid.mat)

    # Remove the lowest tumor stage which will be the reference level for the tumor stage 
    x.train.mat = removeReferenceTumorStage(x.train.mat, levs.rm) 
    x.valid.mat = removeReferenceTumorStage(x.valid.mat, levs.rm) 

    # Remove the reference condition for CN features 
    x.train.mat = dplyr::select(as.data.frame(x.train.mat), -contains(".cn0")) %>% 
                  as.matrix()
                         
    x.valid.mat = dplyr::select(as.data.frame(x.valid.mat), -contains(".cn0")) %>% 
                  as.matrix()
    
    # Finally remove the outcome variables
    x.train.mat = x.train.mat %>% 
           as.data.frame() %>% 
           dplyr::select(!all_of(c(time.var, status.var))) %>% 
           as.matrix()
    
    x.valid.mat = x.valid.mat %>% 
           as.data.frame() %>% 
           dplyr::select(!all_of(c(time.var, status.var))) %>% 
           as.matrix()
    
    ###
    
    # Return object 
    res = list(x.train.mat = x.train.mat,
               y.train = y.train,
               x.valid.mat = x.valid.mat,
               y.valid = y.valid,
               x.train = x.train,
               x.valid = x.valid)
    return(res)
}


#
# Expression data into matrix form
#
expDataToMatrix = function(x){
    x.transp = as.data.frame(t(x))
    header = x.transp[1,]
    x.transp.dat = x.transp[-1,]
    colnames(x.transp.dat) = header
    x.transp.dat.numeric = apply(x.transp.dat, 2, as.numeric)
    rownames(x.transp.dat.numeric) = rownames(x.transp.dat) 
    return(x.transp.dat.numeric)
}

#
# Martrix form back to expression data format 
#

MatrixToExpdata = function(x){
    exp.data = as.data.frame(t(x))
    exp.data = exp.data %>% tibble::rownames_to_column("Participant.ID")
    return(exp.data)
}


#
# Merge data types
#
mergeDataTypes = function(x){
    res.obj = list()
    for (i in 1:length(names(x))){
        end.point = names(x)[i]
        tcga.dataset.splitted.endpoint = tcga.dataset.splitted[[i]] 
        # Merge the training set data across the data types 
        tcga.dataset.splitted.endpoint.train = dplyr::left_join(tcga.dataset.splitted.endpoint$CLIN$train,
                                                                tcga.dataset.splitted.endpoint$EXP$train,
                                                                 by = "Participant.ID")
        tcga.dataset.splitted.endpoint.train = dplyr::left_join(tcga.dataset.splitted.endpoint.train,
                                                                tcga.dataset.splitted.endpoint$CN$train,
                                                                by = "Participant.ID")
        
        # Merge the validation set data across the data types 
        tcga.dataset.splitted.endpoint.valid = dplyr::left_join(tcga.dataset.splitted.endpoint$CLIN$validation,
                                                                tcga.dataset.splitted.endpoint$EXP$validation,
                                                                 by = "Participant.ID")
        tcga.dataset.splitted.endpoint.valid = dplyr::left_join(tcga.dataset.splitted.endpoint.valid,
                                                                tcga.dataset.splitted.endpoint$CN$validation,
                                                                by = "Participant.ID")
        
        res.obj[[end.point]] = list()
        res.obj[[end.point]]$train = tcga.dataset.splitted.endpoint.train
        res.obj[[end.point]]$validation = tcga.dataset.splitted.endpoint.valid
    }
    return(res.obj)
}


