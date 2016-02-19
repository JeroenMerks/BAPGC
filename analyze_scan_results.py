#!/usr/bin/python
"""
Jeroen Merks
Compileerd de resultaten uit de scan in ofwel een Pandas Object of een
csv-bestand. Ook zijn hier zijn alle functies verzameld die te maken hebben
met het berekenen en sorteren van scanresultaten.

"""

import csv
import glob
import os
from itertools import izip

import pandas


def pandalize_csvs(motif_ids, results_folder, merged_results_name):
    print "Compiling scan results..."

    rows = [open(fi).readline().split(",") for fi in
            glob.glob(os.path.join(results_folder + "temp/", "*.csv"))]

    promotor_ids = [line[0] for line in rows]

    # create empty panda dataframe
    dataframe = pandas.DataFrame(index=promotor_ids, columns=motif_ids)
    dataframe.index.name = "promotor_ids"
    dataframe.columns.name = "motif_ids"

    len_motif_ids = len(motif_ids)
    # fill it up row by row
    for row_nr, row in enumerate(rows):
        if len(row[1:]) == len_motif_ids:  # filter out bad results
            dataframe.ix[row_nr] = row[1:]

    # Save as csv for visual inspection
    print "Saving scan results to a comma seperated file..."
    dataframe.to_csv(results_folder + "merged_csv/" + merged_results_name,
                     encoding='utf-8')


def get_best_motif(results_folder, merged_results_name):
    transposed_csv = izip(
        *csv.reader(open(results_folder + merged_results_name, "rU")))
    next(transposed_csv)  # Skip first line with gene id's

    highest_cov_yet = 0
    best_motif = [None]
    for str_row in transposed_csv:
        perc_cov = sum(int(int_row) for int_row in str_row[1:]) / len(
                str_row[1:]) * 100
        if perc_cov > highest_cov_yet:
            highest_cov_yet = perc_cov
            best_motif[0] = str_row[0]

    return best_motif


# Returns list with top x scoring genes
def get_top_100_scoring_genes(results_folder, merged_results_name):
    data = pandas.read_csv(results_folder + merged_results_name,
                           index_col='promotor_ids', dtype='unicode',
                           low_memory=False)

    return data.sum(axis=1).sort_values(ascending=False).head(
        100).keys().tolist()


def get_coregulated_genes(results_folder, results_file, minimum_hits):
    transposed_csv = izip(
        *csv.reader(open(results_folder + results_file, "rU")))
    gene_names = ()
    scores = ()
    coregulated_gene_symbols = []

    for row_nr, row in enumerate(transposed_csv):
        if row_nr == 0:
            gene_names = row
        else:
            scores = row

    for nr, score in enumerate(scores[1:]):
        if int(score) >= minimum_hits:
            coregulated_gene_symbols.append(gene_names[nr + 1].split("_")[-2])

    return coregulated_gene_symbols
