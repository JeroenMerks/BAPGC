from Bio import SeqIO


def main(multiple_fasta_in_file):
    name_longest_gene = None
    length_lonest_gene = 0
    for rec in SeqIO.parse(open(multiple_fasta_in_file, "rU"), "fasta"):
        len_rec = len(str(rec.seq))
        if len_rec > length_lonest_gene:
            length_lonest_gene = len_rec
            name_longest_gene = rec.id.split("_")[-2]

    return name_longest_gene
