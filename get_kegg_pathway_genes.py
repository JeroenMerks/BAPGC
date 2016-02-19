#!/usr/bin/python
"""
Jeroen Merks
Haalt de genen van de opgegeven pathway van KEGG op.

"""
from Bio.KEGG import REST


def main(pathway):
    print "Fetching gene names related to pathway %s from the current KEGG database..." % pathway
    promotor_gene_accessions = []

    pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway

    # iterate through each KEGG pathway file, keeping track of which section
    # of the file we're in, only read the gene in each pathway
    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section

        if current_section == "GENE":
            gene_identifiers, gene_description = line[12:].split("; ")
            gene_id, gene_symbol = gene_identifiers.split()

            if gene_symbol not in promotor_gene_accessions:
                promotor_gene_accessions.append(
                    "promotor_region_" + gene_symbol + "_" + gene_id)

    return promotor_gene_accessions
