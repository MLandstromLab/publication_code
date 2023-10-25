#!/bin/bash

# List files
array=( ACC_workflow.ipynb BLCA_workflow.ipynb BRCA_workflow.ipynb CESC_workflow.ipynb CHOL_workflow.ipynb COAD_workflow.ipynb DLBC_workflow.ipynb ESCA_workflow.ipynb GBM_workflow.ipynb HNSC_workflow.ipynb KICH_workflow.ipynb KIRC_workflow.ipynb KIRP_workflow.ipynb LAML_workflow.ipynb LGG_workflow.ipynb LIHC_workflow.ipynb LUAD_workflow.ipynb LUSC_workflow.ipynb MESO_workflow.ipynb OV_workflow.ipynb PAAD_workflow.ipynb PCPG_workflow.ipynb PRAD_workflow.ipynb READ_workflow.ipynb READ_and_COAD_workflow.ipynb SARC_workflow.ipynb SKCM_workflow.ipynb STAD_workflow.ipynb TGCT_workflow.ipynb THCA_workflow.ipynb THYM_workflow.ipynb UCEC_workflow.ipynb UVM_workflow.ipynb )

for f in ${array[@]};
do 
   jupyter nbconvert --execute --to notebooks --inplace $f
done



