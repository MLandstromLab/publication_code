####################################
## Landstrom Biomarker Project 1  ##
####################################

## imports
import csv
import pathlib
import statistics

## load the data meng
coefficient_fi = "/home/alastairm/Desktop/Landstrom/data/coef.csv"
pvalue_fi = "/home/alastairm/Desktop/Landstrom/data/coef_p.csv"

coefficient_data = []
with open(coefficient_fi, 'r') as infi:
    reader = csv.reader(infi, delimiter=";")
    for row in reader:
        row = [str(x).replace(",",".") for x in row] ## replace european commas
        coefficient_data.append(row)

## load the data meng
pvalue_data = []
with open(pvalue_fi, 'r') as infi:
    reader = csv.reader(infi, delimiter=";")
    for row in reader:
        row = [str(x).replace(",",".") for x in row] ## replace european commas
        pvalue_data.append(row)

## get proteins in lists for index grabbing
column_proteins = coefficient_data[0]
row_proteins = [x[0] for x in coefficient_data]

## Loop over row, find those >0.70 (70%), then grab which column protein the interaction is with
highly_correlated_interactions = []
for row in coefficient_data[1:]:
    for interaction in row[1:]:

        converted_interaction = 0.0
        try: converted_interaction = float(interaction)
        except ValueError: pass ## skip "N/A" in rows

        if converted_interaction >= 0.70:
            column_index = row.index(interaction)
            row_index = coefficient_data.index(row)
            column_protein_interacted = column_proteins[column_index]
            highly_correlated_interactions.append([row[0], converted_interaction, column_protein_interacted, row_index, column_index])

## now get the row/colum indices and grab the P-value for that interaction
for interaction in highly_correlated_interactions:
    row_index = interaction[3]; column_index = interaction[4]
    target_row = pvalue_data[row_index]
    target_value = target_row[column_index]
    interaction.append(target_value)

## form into clean table, save to file
output_rows = [['Protein-1', 'Protein-2', 'Coefficient', 'Pvalue']]
for interaction in highly_correlated_interactions:
    output_rows.append([interaction[0], interaction[2], interaction[1], interaction[5]])
output_fi = "/home/alastairm/Desktop/Landstrom/results/ProteinProteinCoefficients.csv"
with open(output_fi, 'w') as outfi:
    writer = csv.writer(outfi)
    writer.writerows(output_rows)
