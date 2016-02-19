#!/usr/bin/python
"""
Jeroen Merks
Transleer een gegeven multiple-sequence fasta file naar eiwitten.

"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def index_genbank_features(gb_record, feature_type, qualifier):
    answer = dict()
    for index, feature in enumerate(gb_record.features):
        if feature.type == feature_type:
            if qualifier in feature.qualifiers:
                # There should only be one locus_tag per feature, but there
                # are usually several db_xref entries
                for value in feature.qualifiers[qualifier]:
                    if value in answer:
                        pass
                        # print "WARNING - Duplicate key %s for %s features %i and %i" \
                        #       % (value, feature_type, answer[value], index)
                    else:
                        answer[value] = index
    return answer


def locate_cds(top_100_gene_headers, chromosomes_folder, chromosome_file_names):
    cds_seqs = []
    genes_of_interest = [(gene_header, gene_header.split("_")[-1]) for
                         gene_header in top_100_gene_headers]

    for gene_of_interest in genes_of_interest:
        found = False
        gene_header = gene_of_interest[0]
        gene_id = gene_of_interest[1]
        if len(cds_seqs) < 10:

            for chromosome_file_name in chromosome_file_names:
                if not found:
                    gb_records = SeqIO.parse(
                            open(chromosomes_folder + chromosome_file_name,
                                 "r"), "genbank")

                    for gb_record in gb_records:
                        if not found:
                            db_xref_cds_index = index_genbank_features(
                                    gb_record, "CDS", "db_xref")

                            try:
                                gb_feature = gb_record.features[
                                    db_xref_cds_index["GeneID:" + gene_id]]
                                cds_seqs.append((gene_header,
                                                 gb_feature.extract(
                                                         gb_record.seq)))  # (gene_header, cds_seq)
                                print "Found CDS site of " + gene_header
                                found = True
                                break
                            except KeyError:  # sequence not in chromosome
                                pass
                        else:
                            break
                else:
                    break
            if not found:
                print gene_header + " has no known CDS sites..."

    return cds_seqs


def make_prot_rec(gene_header, cds_seq):
    """Maakt een BioPython SeqRecord van de getransleerde DNA-sequentie."""
    # slice off stop codon *
    return SeqRecord(seq=cds_seq.translate()[:-1], id=gene_header,
                     description="translation of CDS of gene, default table")


def translate_cds_seqs(cds_seqs):
    proteins = []
    for gene_header, cds_seq in cds_seqs:
        proteins.append(make_prot_rec(gene_header, cds_seq))

    return proteins


def extract_gene_headers(multiple_fasta_dna_in_file):
    gene_headers = []
    handle = open(multiple_fasta_dna_in_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        gene_headers.append(record.id)

    handle.close()
    return gene_headers


def main(multiple_fasta_dna_in_file, multiple_fasta_prot_out_file,
         chromosomes_folder, chromosome_file_names):
    in_handler = open(multiple_fasta_dna_in_file, "rU")
    out_handler = open(multiple_fasta_prot_out_file, "w+")

    gene_headers = extract_gene_headers(multiple_fasta_dna_in_file)

    print "Collecting CDS sites of genes..."
    cds = locate_cds(gene_headers, chromosomes_folder, chromosome_file_names)

    print "Translating CDS sites to proper proteins..."
    proteins = translate_cds_seqs(cds)

    for protein in proteins:
        SeqIO.write(protein, out_handler, "fasta")

    in_handler.close()
    out_handler.close()
