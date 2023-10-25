library(dplyr)
source("../biomartQuery.R")

# 
# Function for obtaining Copy number data from TCGA stored
# on server. Takes cancer type (following TCGA acronyms) as input
#
getTCGACopyNumberData = function(cancer.type){
  
  # Cancer type with prefix
  ctype = paste0("TCGA-", cancer.type)
  
  # Upstream path
  tcga.cn.path.up = file.path("/home/data/TCGA/Copy_Number_Variation/", ctype)
  
  # Find path
  path.cn = list.files(path = tcga.cn.path.up, 
                       pattern = "GeneLevelCopyNumberScores.RData", 
                       recursive = T) 
  
  # Full path 
  full.path.cn = file.path(tcga.cn.path.up, path.cn )
  
  # Cn data 
  cn.data = load(full.path.cn)
  cn.data = get(cn.data)
  
  # Create the data object
  data.obj = list("Clinical.data" = as.data.frame(colData(cn.data)),
                  "Data" =  assay(cn.data),
                  "Annotation" = as.data.frame(rowData(cn.data)))
  
  return(data.obj)                          
}


#
# Add minimal annotation to gene expression data 
# We will add the gene symbols 

addAnnotation = function(ensembl.ids, annotation.file) {

    # Read in the annotation file 
    id.to.symbol = read.csv(annotation.file)

    # Do left join
    annotation.data = as.data.frame(ensembl.ids)
    colnames(annotation.data) = "gene_id"
    annotation.data = dplyr::left_join(annotation.data, id.to.symbol, by = "gene_id")

    annotation.data$gene_name[annotation.data$gene_name == ""] = NA
    
    # Few ENS ids seem to have multiple matching symbols. We will compress
    # these with ";"
    split.data = split(annotation.data, annotation.data$gene_id)
    gene.symbols.compressed = unlist(map(split.data, function(x){paste(x$gene_name, collapse = ";")}))

    annotation.data.comp = data.frame("Gene.Symbol" = gene.symbols.compressed,
                                    "Gene.ID" = names(split.data))
  
    # Make sure the order is correct 
    annotation.data.comp = annotation.data.comp[match(ensembl.ids, annotation.data.comp$Gene.ID),]

    return(annotation.data.comp)
}



# 
# Function for obtaining gene expression from TCGA stored
# on server. Takes cancer type (following TCGA acronyms) as input
# The data is non-normalised count data 
getTCGAExpressionData = function(cancer.type, annotation.file){
  
  # Cancer type with prefix
  ctype = paste0("TCGA-", cancer.type)
  
  # Upstream path
  tcga.exp.path.up = file.path("/home/data/TCGA/Gene_Expression_Quantification/", ctype)
  
  # Find path
  path.exp = list.files(path = tcga.exp.path.up, 
                        pattern = "HTSeq-Counts.RData", 
                        recursive = T)   
  
  # Full path 
  full.path.exp = file.path(tcga.exp.path.up, path.exp)
  
  # Exp data 
  exp.data = load(full.path.exp)
  exp.data = get(exp.data)

  
  # Return an object with count data and clinical data
  data.obj = list("Clinical.data" = as.data.frame(colData(exp.data)),
                  "Data" = assay(exp.data))
  
  # Add minimal annotation data 
  data.obj$Annotation = addAnnotation(rownames(data.obj$Data), annotation.file)
  
  return(data.obj)
  
}

# 
# Function of obtaining masked mutation data from TCGA stored on the server 
#
getTCGAMutationData = function(cancer.type){

  # Cancer type with prefix
  ctype.pattern = paste0("TCGA.", cancer.type, ".mutect")
  
  # Upstream path
  tcga.mut.path.up = file.path("/home/data/TCGA/Masked_Somatic_Mutation")
  
  # Find path
  path.mut = list.files(path = tcga.mut.path.up, 
                       pattern = ctype.pattern,
                       full.name = T)

  # Read file 
  mut.data = read.csv(path.mut)

  return(mut.data)
}


#
# Function for obtaining cancer specific clinical data
#
#
getClinData = function(cancer.type){
  
  # Path to clinical data
  path.to.clin = "/home/data/TCGA/Clinical_Data"
  
  # File name 
  fname = paste0("TCGA-", cancer.type, "_clinical.csv")
  
  # Read in data 
  clin.data = read.csv(file.path(path.to.clin, fname))
  
  return(clin.data)
}


#
# Function of obtaining clinical end point data 
# 
#
getClinEndpointData = function(cancer.type){
  
  # Read in end point data 
  clin.endpoints = read.csv("/home/data/TCGA/TCGA_CDR_clinical_endpoints.csv") 
  
  # Restrict analysis to cancer type of interest 
  clin.endpoints.can = dplyr::filter(clin.endpoints, type == cancer.type)
  
  # Return 
  return(clin.endpoints.can)
  
}

