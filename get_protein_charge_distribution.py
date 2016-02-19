#!/usr/bin/python
"""
Christian van der Niet
21-01-2016
programma berekent percentage van hydrofobe, hydrofiele en neutrale amino zuren
en geeft de verhousing weer welke amino zuren in een seqentie zitten in een staafdiagram
"""
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO


def main(protein_multiple_fasta, bar_plot_pic):
    # lijst met alle amino zuren
    letter = ["A", "R", "N", "D", "C", "F", "Q", "E", "G", "H", "I", "L", "K",
              "M", "P", "S", "T", "W", "Y", "V"]
    # lijst met alle ladingen
    lading = ["0", "+", "0", "-", "0", "0", "0", "-", "0", "+", "-", "-", "+",
              "0", "0", "0", "0", "0", "0", "0"]

    lijstpecentages = []
    lijstexactehoeveelheid = []
    lijstnamen = []

    for record in SeqIO.parse(protein_multiple_fasta, "fasta"):

        # naam van eiwit
        name = record.id.rstrip()
        # lijst met alle hoeveelheden
        hoeveelheid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0]
        # totaal aantal, positief, neutraal, negatief geladen
        lading2 = [0, 0, 0, 0]
        lijstnamen.append(
            "%s positief %s neutraal %s negatief" % (name, name, name))

        for aa in list(record.seq):

            # check lading met if statement
            if lading[letter.index(aa)] == "+":
                lading2[0] += 1
                lading2[1] += 1
            elif lading[letter.index(aa)] == "0":
                lading2[0] += 1
                lading2[2] += 1
            elif lading[letter.index(aa)] == "-":
                lading2[0] += 1
                lading2[3] += 1

            # check hoeveelheid met if statement
            hoeveelheid[letter.index(aa)] += 1

        # lading:
        # percentage waardes van de positief neutraal en negatief
        positief = float(lading2[1]) * (100 / float(float(lading2[0])))
        neutraal = float(lading2[2]) * (100 / float(float(lading2[0])))
        negatief = float(lading2[3]) * (100 / float(float(lading2[0])))
        print "percentage = %s positief %s neutraal %s negatief" % (
        positief, neutraal, negatief)

        # zet in lijst
        lijstpecentages.append("%.1f%s" % (positief, "%"))
        lijstpecentages.append("%.1f%s" % (neutraal, "%"))
        lijstpecentages.append("%.1f%s" % (negatief, "%"))
        lijstexactehoeveelheid.append(lading2[1])
        lijstexactehoeveelheid.append(lading2[2])
        lijstexactehoeveelheid.append(lading2[3])

        # hoeveelheid:
        # get_bar_plot(hoeveelheid, letter, bar_plot_name)

    plot_verhouding(lijstnamen, lijstpecentages, lijstexactehoeveelheid,
                    bar_plot_pic)


def get_bar_plot(hoeveelheid, letter, bar_plot_pic):
    # geef variabelen op
    n = 20
    X = np.arange(n)
    Y = hoeveelheid
    my_xticks = letter
    # maak coole x as met letters
    plt.xticks(X, my_xticks)
    # maak een bar plot
    plt.bar(X, Y)
    # maak mooie titeltjes en label de assen
    plt.title('Aminozuursamenstelling')
    plt.xlabel('Aminozuur')
    plt.ylabel('Frequentie')

    # sla de afbeelding op
    plt.savefig(bar_plot_pic)


def plot_verhouding(lijstnamen, lijstpecentages, lijstexactehoeveelheid,
                    plot_file_name):
    plt.clf()
    # print lijstnamen
    # print lijstpecentages
    # print lijstexactehoeveelheid
    # geef variabelen op
    n = len(lijstexactehoeveelheid)
    # print len(lijstexactehoeveelheid)
    X = np.arange(n)
    # print X
    Y = lijstexactehoeveelheid
    # print Y
    my_xticks = lijstpecentages
    # print my_xticks
    # maak coole x as met letters
    plt.xticks(X, my_xticks)
    # maak een bar plot
    plt.bar(X, Y)
    # maak mooie titeltjes en label de assen
    plt.title('aminozuursamenstelling')
    plt.xlabel("     ".join(lijstnamen))
    plt.ylabel('hoeveelheid')
    # verander font
    plt.rcParams.update({'font.size': 6})
    # sla de afbeelding op
    plt.savefig(plot_file_name)
