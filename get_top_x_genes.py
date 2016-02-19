#!/usr/bin/python
"""
Jeroen Merks
Leest een grote multiple-sequence fasta bestand in en slaat de eerste top_x aantal records op in een apparte fasta.

"""
from Bio import SeqIO


def main(multiple_fasta_in_file, top_genes_folder, top_x):
    out_file = top_genes_folder + "top_%i_genes.fasta" % top_x
    output_handle = open(out_file, "w+")
    recs = []
    for nr, rec in enumerate(
            SeqIO.parse(open(multiple_fasta_in_file, "rU"), "fasta")):
        recs.append(rec)
        if nr == top_x:
            SeqIO.write(recs, output_handle, "fasta")
            output_handle.close()

    return out_file
