"""
Christian van der Niet
programma maakt een msa uit gegeven input faste een msa


"""
from Bio import AlignIO


def main(alignment_in):
    # print align file
    align = AlignIO.read(alignment_in, "clustal")
    # geconserverde nucleotides
    overeenkomend = (align.format("clustal").count("*"))
    # lengte langste sequentie
    first_gene = (align.format("clustal").split()[5])
    # hoeveelheid keer een seqentie in een nieuwe regel staat
    nr_sequences = (align.format("clustal").split()).count(first_gene)
    alignment = (align.format("clustal").split())
    # loop extracts the first sequence from the msa
    seqence = ""

    for x in range(nr_sequences):
        pos = alignment.index(first_gene) + 1
        seqence += (alignment[pos])
        alignment.pop(alignment.index(first_gene))

    # bereken percentage
    percentage_conserved = float(
        float(100) / float(len(seqence)) * float(overeenkomend))

    return percentage_conserved
