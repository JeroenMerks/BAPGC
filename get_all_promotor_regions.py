#!/usr/bin/python
"""
Jeroen Merks
Haal alle promotorregionen van alle bekende "gene" accesions in de GenBank
bestanden van het genoom van de mens op. De promotorregio begint 6000
basenparen voor de start van het gen en 500 na.

"""

from Bio import SeqIO


def main(chromosomes_folder, chromosome_names, promotors_folder):
    for chromosome_file_name in chromosome_names:

        # hs_ref_GRCh38.p2_chr1.gbk -> chr1
        chromosome_nr = chromosome_file_name.split("_")[-1][:-4]

        print "Reading chromosome %s..." % chromosome_nr
        promotors_output_loc = promotors_folder + "promotors_" + chromosome_nr + ".fasta"
        promotor_seqs_fasta = open(promotors_output_loc, "w+")

        print "Extracting and saving promotor regions of all known gene sites" \
              " of chromosome '%s'..." % chromosome_nr
        for genome in SeqIO.parse(chromosomes_folder + chromosome_file_name,
                                  "genbank"):
            for feature in genome.features:
                if feature.type == "gene":
                    gene_name = feature.qualifiers["gene"][0]
                    gene_id = feature.qualifiers["db_xref"][0].split(":")[1]
                    prom_reg_start = feature.location.start - 6000
                    prom_reg_end = feature.location.start + 500

                    if prom_reg_start > 0:
                        promotor_seqs_fasta.write(
                                ">promotor_region_" + gene_name + "_" + gene_id + "\n" + str(
                                        genome.seq[prom_reg_start:
                                        prom_reg_end]) + "\n")

        promotor_seqs_fasta.close()
