#!/usr/bin/python
"""
Jeroen Merks
Haalt de genen op van promotorregionen op het genoom die hoog gescoort hebben.

"""
from Bio import SeqIO


def main(chromosomes_folder, chromosome_file_names, promotor_ids,
         specific_genes_fasta_name):
    # Extract gene id's from the promotor fasta id's
    # [">promotor_region_HMGCL_3155", ">promotor_region_OXCT2_64064"] -> ["3155", "64064"]
    gene_ids = [prom_id.split("_")[-1] for prom_id in promotor_ids]

    specific_genes_fasta_file = open(specific_genes_fasta_name, "w+")

    for chromosome_file_name in chromosome_file_names:

        for genome in SeqIO.parse(chromosomes_folder + chromosome_file_name,
                                  "genbank"):
            for feature in genome.features:
                if feature.type == "gene":

                    gene_id = feature.qualifiers["db_xref"][0].split(":")[1]
                    if gene_id in gene_ids:
                        gene_ids.remove(gene_id)

                        gene_name = feature.qualifiers["gene"][0]
                        gene_reg_start = feature.location.start
                        gene_reg_end = feature.location.end

                        specific_genes_fasta_file.write(
                            ">gene_" + gene_name + "_" + gene_id + "\n" + str(
                                    genome.seq[
                                    gene_reg_start:gene_reg_end]) + "\n")

    specific_genes_fasta_file.close()
