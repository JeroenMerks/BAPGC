#!/usr/bin/python
"""
Jeroen Merks
Haalt specifieke promotorregionen op die gerelateerd zijn aan de opgegeven
pathway.

"""
import os

from Bio import SeqIO


def main(pathway_gene_names, promotor_regions_file, all_promotors_folder):
    output_file = open(promotor_regions_file, "w+")

    promotor_file_names = [filename for filename in
                           os.listdir(all_promotors_folder) if
                           filename[-6:] == '.fasta']
    for promotor_file in promotor_file_names:
        handle = all_promotors_folder + promotor_file

        for promotor_region_gene in SeqIO.parse(open(handle, "rU"), "fasta"):
            for pathway_gene_name in pathway_gene_names:
                if pathway_gene_name == promotor_region_gene.id.rstrip():
                    pathway_gene_names.remove(pathway_gene_name)
                    SeqIO.write(promotor_region_gene, output_file, "fasta")
