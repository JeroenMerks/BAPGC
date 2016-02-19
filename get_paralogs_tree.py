#!/usr/bin/python
"""
 Rick de Graaf
 s1072043
 getParalogs.py

"""

# laden van modulen
import re

import mygene
import requests
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from ete3 import Tree


# Inlezen en extraheren van de symbolen uit het top_50_genes.fasta bestand
def extractGenesFasta(top_50_genes_file_name):
    # openen van bestand
    handler = open(top_50_genes_file_name, "rU")
    genen_list = []
    # extraheren van symbols
    for record in SeqIO.parse(handler, "fasta"):
        rec = record.id
        gene_symbol = rec.split("_")[-2]
        genen_list.append(gene_symbol)
    handler.close()
    # filteren van onbruikere genen
    filtered_list = [rec for rec in genen_list if "LOC" not in rec]
    gene_symbols = [rec for rec in filtered_list if "orf" not in rec]
    # teruggeven van symbolen
    return gene_symbols


# Selecteren van de 10 beste genen deze worden vervolgens geconverteerd van
# symbols naar ensemble id's
def getEnsembleNames(gene_list):
    mg = mygene.MyGeneInfo()
    # selecteren van de 10 beste genen
    genes_selection = gene_list[0:10]
    # converten van symbols naar ensemble id' doormiddel van de mygene module
    convert_genes = mg.querymany(genes_selection, scopes='symbol',
                                 fields='ensembl.gene', species='human')
    genes = []
    # controleer dat er geen lege velden zijn
    for item in convert_genes:
        variable = item.get('ensembl.gene', None)
        if variable is None:
            pass
            # print("Empty field, skip")
        else:
            genes.append(variable)
    # teruggeven van genen
    return genes


# Verkrijven paralogen van de genen doormiddel van de Ensemble rest api
# Na het verkrijgen van de paralogen wordt de nuttige informatie geselecteerd
# en weggeschreven naar een FASTA bestand
def getParalogs(gene_symbols, paralogs_file):
    # Maken van FASTA bestand
    handler = open(paralogs_file, "w+")
    print("Finding paralogs of top 50 genes...")
    # Loopen door de lijst met genen
    for gene_symbol in gene_symbols:
        # toegang verkrijgen tot de Ensemble rest API
        server = "http://rest.ensembl.org"
        ext = "/homology/symbol/human/%s?type=paralogues;sequence=cdna;target_species=homo_sapiens;cigar_line=0;aligned=0" % (
            gene_symbol)
        r = requests.get(server + ext, headers={"Content-Type": "text/xml"})
        paralog_info = r.text
        # Alleen de paralogen selecteren
        target = re.findall('<target id=(.*)\/>', paralog_info)
        # Selecteren van de best scorende paralogen
        scorelist = []
        # Extraheren van de scores
        for item in target:
            score = re.findall('perc_id="(.*)\" perc_pos', item)
            scorelist.append(score[0])
        scorelist.sort(reverse=True)
        scorelist = scorelist[:4]
        # Extraheren van de informatie om het FASTA bestand te bouwen
        # Id van het gen, sequentie en score
        for item in target:
            score = re.findall('perc_id="(.*)\" perc_pos', item)
            if score[0] in scorelist:
                gene_symbol = re.findall('"(.*)\" perc_id', item)
                seq = re.findall('seq="(.*)\" species=', item)
                # Weggschrijven van de informatie naar het FASTA bestand
                fasta = ">%s_%s\n%s\n" % (gene_symbol[0], score[0], seq[0])
                handler.write(fasta)
    handler.close()


# Fasta bestand wordt ingeladen, hier wordt vervolgens een ClustalW2 alignment
# op uitgevoerd. Vervolgens wordt er een boom gecreerd die wordt
# weggschreven naar een PNG bestand
def builtTree(phylo_tree_pic, paralogs_file):
    print("Aligning top 50 genes for phylogenetic tree...")
    # maken van alignment
    clustalw_cline = ClustalwCommandline("clustalw", infile=paralogs_file)
    stdout, stderr = clustalw_cline()
    # importeren van boom bestand
    tree = Tree(paralogs_file[:-6] + ".dnd")
    # bouwen en weggschrijven van boom
    tree.render(phylo_tree_pic)


# main functie voor het aanroepen van de overige functies        
def main(top_50_genes_file_name, paralogs_file, phylo_tree_pic):  # argv
    # instellen van bestandsnaam
    # verkrijgen van gen symbolen
    gene_symbols_list = extractGenesFasta(top_50_genes_file_name)
    # verkrijgen van Ensemble namen
    ensembl_gene_names = getEnsembleNames(gene_symbols_list)
    # verkrijgen van paralogen
    getParalogs(ensembl_gene_names, paralogs_file)
    # instellen van png bestandnaam
    # bouwen van boom
    builtTree(phylo_tree_pic, paralogs_file)

# Extra informatie:
# =======================
#
# Dit script leest een FASTA bestand in met de beste meegereguleerde genen, en
#  pakt
# vervolgens de 10 beste genen, deze genen worden vervolgens omgezet naar
# ENSEMBLE id's.
# Van deze id's worden vervolgens de bekende paralogen genen opgehaald uit de
# ENSEMBLE database.
# Hiervan worden vervolgens de vier beste genen gebruikt om eerst een alignment
#  te maken
# en vervolgens om een boom te bouwen die wordt weggeschreven als PNG bestand.
