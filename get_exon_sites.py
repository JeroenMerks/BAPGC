#!/usr/bin/python
"""
Koen van Diemen
Haal de locaties van de intronen van de gegeven genen op en bereken op basis daarvan de verhoudingen van intron- en exon
lengte. Deze verhoudingen worden in een scatterplot opgeslagen.

"""
import matplotlib.pyplot as plt
from Bio import SeqIO


def Get_GeneID_and_Length(top_20_genes_file):
    records = list(SeqIO.parse(top_20_genes_file, "fasta"))
    GeneIDS = []

    for record in records:
        gene_id = record.id.split("_")[-1]
        seq = len(record.seq)

        # Alle gen namen in een lijst, met daarin de lenghte van het gen. Lengthe wordt
        # Verder op gebruikt om de intron lengte te berekenen.
        GeneIDS.append("GeneID:" + gene_id)
        GeneIDS.append(seq)

    return GeneIDS


def Get_Exon_intron_length(chromosomes_folder, chromosome_names, GeneIDS):
    intronen = []
    exonen = []

    for chromosome_file_name in chromosome_names:

        # print "Finding intron sites of the genes on chromosome %s..." % chromosome_nr

        for genome in SeqIO.parse(chromosomes_folder + chromosome_file_name,
                                  "genbank"):
            for feature in genome.features:

                # Als het gen de feature "exon" bevat in de genbank en dit is 1 van de top 50
                # genen is er dus een exon gevonden in het gen, hiervan wordt de lengte bepaald.
                if feature.type == "exon" and feature.qualifiers['db_xref'][
                    0] in GeneIDS:
                    # gene_name = feature.qualifiers["gene"][0]

                    # Gene ID en genaam worden nog een keer weergegeven.
                    # print(feature.qualifiers['db_xref'][0], "Gene_name:" + gene_name)

                    # lengte van het gen wordt opgehaald.
                    gene_index = GeneIDS.index(feature.qualifiers['db_xref'][0])
                    gene_length = GeneIDS[gene_index + 1]

                    # Start en stop positie van het exon wordt geselecteerd.
                    start_position = feature.location.start
                    end_position = feature.location.end

                    # Cumalatieve intron en exon lengte worden bepaald voor een bepaald gen.
                    exon_length = end_position - start_position
                    intron_length = gene_length - exon_length
                    # print("exon_length:" + str(exon_length), "Sequence_length:" + str(gene_length),
                    #       "intron_length:" + str(intron_length))

                    intronen.append(intron_length)
                    exonen.append(exon_length)

    return intronen, exonen


def get_scatter_plot(intronen, exonen, scatter_plot_pic_name):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(exonen, intronen)

    # zet de namen van x en y as
    ax1.set_xlabel('Exon lengte')
    ax1.set_ylabel('Intron lengte')

    # als nodig is zet deze assen voor betere plaatjes.
    ax1.set_ylim([-2, max(intronen) + 10])
    ax1.set_xlim([0, max(exonen) + 10])

    # save scatterplot to file
    # plt.show()
    print "Generating scatter plot..."
    fig.savefig(scatter_plot_pic_name)


def main(chromosomes_folder, chromosome_file_names, top_20_genes_file,
         plot_pic_name):
    genes = Get_GeneID_and_Length(top_20_genes_file)
    Gen_length = Get_Exon_intron_length(chromosomes_folder,
                                        chromosome_file_names, genes)
    get_scatter_plot(Gen_length[0], Gen_length[1], plot_pic_name)
