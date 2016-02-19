#!/usr/bin/python
"""
Jeroen Merks
Haalt het plaatje van de pathway op van KEGG.

"""
import urllib


def main(met_pathway_picture_folder, met_pathway_code):
    met_pathway_picture_file = met_pathway_picture_folder + met_pathway_code + ".png"
    handler = open(met_pathway_picture_file, 'w+')
    handler.write(urllib.urlopen(
        "http://www.genome.jp/kegg/pathway/hsa/%s.png" % met_pathway_code).read())
    handler.close()

    return met_pathway_picture_file
