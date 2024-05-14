####################################
## Landstrom Biomarker Project 1  ##
####################################

## imports
import csv
import pathlib
import statistics

## load the data meng
disease_fi = "/home/alastairm/Desktop/Landstrom/results/diseaseValues.csv"
control_fi = "/home/alastairm/Desktop/Landstrom/results/controlValues.csv"

disease_data = []
with open(disease_fi, 'r') as infi:
    reader = csv.reader(infi)
    for row in reader:
        toAppend = row[:-14] ## remove/ignore the "phenotype" data at the end of each row
        disease_data.append(toAppend)

control_data = []
with open(control_fi , 'r') as infi:
    reader = csv.reader(infi)
    for row in reader:
        toAppend = row[:-14] ## remove/ignore the "phenotype" data at the end of each row
        control_data.append(toAppend)

top70_proteins = ["ESM.1","MK","SYND1","FGF.BP1","IL6","TFPI.2","HGF","ANXA1","WFDC2","MetAP.2","TNFRSF6B","RET","DLL1","CAIX","CXCL13","TNFRSF4","WISP.1","FR.alpha","CD70","EGF","RSPO3","VEGFA","TXLNA","LYN","S100A11","PVRL4","hK14","ADAM.8","VIM","AREG","GZMB","KLK13","hK11","VEGFR.3","CD27","FADD","SCAMP3","CD207","LY9","ABL1","IFN.gamma.R1","MSLN","CTSV","FURIN","MUC.16","CDKN1A","CD48","FASLG","CEACAM1","ICOSLG","CD160","ERBB3","ITGB5","FR.gamma","Gal.1","CPE","TCL1A","S100A4","ERBB4","hK8","FCRLB","PPY","TGFR.2","IGF1R","TRAIL","MAD.homolog.5","GPC1","MIC.A.B","WIF.1","GPNMB"]

master_data = []
for protein in top70_proteins:                                              ## For every protein we want to look at
    for target in disease_data[0]:                                          ## Loop over the protein column names in disease data
        if protein == target:                                               ## Match between current top70 protein and column name protein
            columnIndex = disease_data[0].index(target)                     ## Get column index for the matched top70/column name protein
            for row in disease_data[1:]:                                    ## Loop over every row (sample/individual) in disease data (skipping protein column name/first row)
                master_data.append(["Disease", protein, row[columnIndex]])  ## Add the NPX reading at the protein-index, from each row (sample/individual)

for protein in top70_proteins:                                              ## Repeat the same procedure for control_data 
    for target in control_data[0]:
        if protein == target:
            columnIndex = control_data[0].index(target)
            for row in control_data[1:]:
                master_data.append(["Control", protein, row[columnIndex]])

## Write master_data to file, for plotting in R
output_fi = "/home/alastairm/Desktop/Landstrom/results/boxplotData.csv"
with open(output_fi, 'w') as outfi:
    writer = csv.writer(outfi)
    writer.writerows(master_data)
