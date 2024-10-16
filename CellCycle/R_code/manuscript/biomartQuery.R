library("biomaRt")

getGeneSymbols = function(ens.ids) {
  
  # 
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Get genes names                        
  res = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                                     filters = "ensembl_gene_id", 
                                     values = ens.ids, 
                                     mart = ensembl)
  
  return(res)
}


