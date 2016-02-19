#!/usr/bin/python
"""
Jeroen Merks
De pipeline van het analyseprocess.

"""
import multiprocessing
import os
import pickle
import sys

from Bio import SeqIO

import analyze_scan_results
import get_all_promotor_regions
import get_exon_sites
import get_kegg_pathway_genes
import get_met_pathway_picture
import get_motifs_JASPAR
import get_msa
import get_name_longest_gene
import get_paralogs_tree
import get_percentage_conserved
import get_protein_charge_distribution
import get_specific_genes
import get_specific_promotor_regions
import get_top_x_genes
import motif_scanner
import pdf_generator
import translate_dna

# Installeer dependencies, getest op Ubuntu 15.10
os.system("bash ./install_dependencies.sh")

# Alle relatieve locaties van bestanden en mappen binnen de mappenstructuur
results_folder = "data/results/"
chromosomes_folder = "data/genome/GRCh38/"
all_motifs_folder = "data/motifs/motifs_JASPAR/"
all_promotors_folder = "data/promotors/whole_genome/"
pathway_promotors_folder = "data/promotors/pathway_specific/"
pathway_scan_results_folder = results_folder + "pathway_specific/"
whole_genome_results_folder = results_folder + "whole_genome/"
msa_file = results_folder + "msa/msa.aln"
phylo_tree_folder = results_folder + "phylo_tree/"
paralogs_file = phylo_tree_folder + "paralogs.fasta"
phylo_tree_pic = phylo_tree_folder + "phylogenetic_tree_paralogs.png"
intron_exon_scatter_plot_pic = results_folder + "intron_exon_scatter_plot/intron_exon_scatter_plot.png"
top_100_genes_folder = results_folder + "genes/"
top_100_genes_file = top_100_genes_folder + "top_100_genes.fasta"
proteins_folder = results_folder + "proteins/"
proteins_file = proteins_folder + "top_10_proteins.fasta"
proteins_charge_distribution_pic = proteins_folder + "protein_charge_distribution.png"
pathway_promotors_file = pathway_promotors_folder + "pathway_promotors.fasta"
met_pathway_picture_folder = results_folder + "met_pathway_picture/"
pickle_file = results_folder + "pickle/pickle.p"

# Alle JASPAR motifs in .pfm formaat.
all_motif_ids = [filename for filename in os.listdir(all_motifs_folder) if
                 filename[-4:] == '.pfm']
# Alle humane chromosomen in GenBank .gbk formaat.
chromosome_file_names = [filename for filename in os.listdir(chromosomes_folder)
                         if filename[-4:] == '.gbk']

CPU_count = multiprocessing.cpu_count() - 1  # Stel het aantal processorkernen
# vast en gebruik deze waar mogelijk.
# Houdt 1 kern over voor de OS

# De parameters m.b.t. het scannen van motifs.
pseudocount = 0.00001
pvalue = 0.0001

studentnummer = "s1072614"  # Studentnummer van onderzoeker, pas aan
# indien nodig.
auteurs = "Melissa van Wieringen, Jeroen Merks, Koen van Diemen, Rick " \
          "de Graaf en Christian van der Niet"

# Alle namen van motifs die volgens Jaspar CORE gerelateerd zijn aan de mens.
# MA0528.1 en MA0050.2 zijn weggelaten, omdat deze
# (ook na de normalisatiestappen) significant vaker hitten dan de andere motifs
all_human_motif_ids = ["MA0017.1", "MA0025.1", "MA0028.1", "MA0030.1",
                       "MA0031.1", "MA0032.1", "MA0033.1", "MA0042.1",
                       "MA0043.1", "MA0048.1", "MA0051.1", "MA0056.1",
                       "MA0057.1", "MA0059.1", "MA0066.1", "MA0069.1",
                       "MA0070.1", "MA0071.1", "MA0072.1", "MA0073.1",
                       "MA0074.1", "MA0077.1", "MA0081.1", "MA0084.1",
                       "MA0090.1", "MA0091.1", "MA0101.1", "MA0107.1",
                       "MA0115.1", "MA0119.1", "MA0124.1", "MA0130.1",
                       "MA0131.1", "MA0139.1", "MA0149.1", "MA0138.2",
                       "MA0112.2", "MA0155.1", "MA0159.1", "MA0161.1",
                       "MA0163.1", "MA0462.1", "MA0465.1", "MA0466.1",
                       "MA0468.1", "MA0470.1", "MA0471.1", "MA0473.1",
                       "MA0475.1", "MA0476.1", "MA0477.1", "MA0478.1",
                       "MA0479.1", "MA0481.1", "MA0484.1", "MA0486.1",
                       "MA0488.1", "MA0489.1", "MA0490.1", "MA0491.1",
                       "MA0492.1", "MA0495.1", "MA0496.1", "MA0497.1",
                       "MA0501.1", "MA0502.1", "MA0504.1", "MA0506.1",
                       "MA0507.1", "MA0508.1", "MA0510.1", "MA0511.1",
                       "MA0513.1", "MA0516.1", "MA0517.1", "MA0523.1",
                       "MA0524.1", "MA0525.1", "MA0526.1", "MA0527.1",
                       # "MA0528.1",
                       "MA0007.2", "MA0102.3", "MA0024.2", "MA0154.2",
                       "MA0162.2", "MA0076.2", "MA0258.2", "MA0148.3",
                       "MA0036.2", "MA0037.2", "MA0114.2",  # "MA0050.2",
                       "MA0058.2", "MA0052.2", "MA0105.3", "MA0060.2",
                       "MA0014.2", "MA0079.3", "MA0083.2", "MA0137.3",
                       "MA0144.2", "MA0140.2", "MA0003.2", "MA0106.2",
                       "MA0093.2", "MA0095.2", "MA0103.2", "MA0592.1",
                       "MA0593.1", "MA0595.1", "MA0596.1", "MA0597.1",
                       "MA0598.1", "MA0599.1", "MA0600.1"]


# Indien de gebruiker de help-vlag meegeeft, geef de gebruikersinformatie weer.
def show_usage_information():
    print()
    print("""Dit is een pipeline voor het ontdekken van genen buiten een
    opgegeven pathway waarvan sterke aanwijzingen bestaan dat deze door
    dezelfde transcriptiefactoren worden gereguleerd. De gevonden genen
    worden op verschillende manieren verder
 onderzocht.""")
    print()
    print("""Gebruik:
    python(2) start_pipeline.py '[PATHWAY_CODE]'

    PATHWAY_CODE, de KEGG pathway code van de metabole route

    Voorbeeld:
    python start_pipeline.py 'hsa04150'""")
    print()
    sys.exit()


# Functie die folders recursief legen.
def delete_previous_results(results_folder):
    empty_dirs = []

    def delete_files(dir_list, dir_path):
        for f in dir_list:
            if f is not ".keep":
                os.remove(dir_path + "/" + f)

    def clean_directory(dir_entry):
        delete_files(dir_entry[2], dir_entry[0])
        empty_dirs.insert(0, dir_entry[0])

    tree = os.walk(results_folder)
    for directory in tree:
        clean_directory(directory)


# De main functie
def main(argv):
    if len(argv) >= 2:
        if argv[1] == "-h" or argv[1] == "--h" or argv[1] == "-help" or argv[
            1] == "--help":
            show_usage_information()

    met_pathway_code = argv[1]

    delete_previous_results(results_folder)

    # Download de assembly van het menselijke genoom van NCBI indien de folder
    # van het chromosoom leeg is.
    if len(os.listdir(chromosomes_folder)) == 0:
        retval = os.getcwd()  # Sla de working directory op, zodat hier later
        #  weer naar terug kan
        os.chdir(chromosomes_folder)  # cd naar de folder waar het chromosoom
        #  terecht moet komen.
        os.system('wget -r -nd -X /genomes/H_sapiens/ARCHIVE,'
                  '/genomes/H_sapiens/ARCHIVE '
                  '-t 45 --accept="hs_ref_GRCh38.p2_chr*.gbk.gz" '
                  '"ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens*"')
        print "Extracting genome files..."
        os.system('gunzip *.gz')  # Pak de GZ-archieven van de chromosomen uit
        os.chdir(retval)  # cd terug naar de working directory.

    # print "Downloading pathway picture %s from KEGG..." % met_pathway_code
    met_pathway_pic = get_met_pathway_picture.main(met_pathway_picture_folder,
                                                   met_pathway_code)

    # Download alle motifs van JASPAR CORE
    get_motifs_JASPAR.main(all_human_motif_ids, all_motifs_folder)

    # Vraag bij KEGG op welke genen betrokken zijn bij de opgegeven pathway.
    pathway_related_promotor_ids = get_kegg_pathway_genes.main(met_pathway_code)

    # Indien nog niet alle promotorregionen gedownload zijn
    if len(os.listdir(all_promotors_folder)) == 0:
        # Haal deze op uit de GenBank-files van het genoom.
        get_all_promotor_regions.main(chromosomes_folder, chromosome_file_names,
                                      all_promotors_folder)

    # Alle bestanden met een .fasta-extensie in de map voor promotors, worden
    # gezien als multiple-fasta bestanden van
    # alle promotorreionen van de mens.
    all_promotors_file_names = [fn for fn in os.listdir(all_promotors_folder) if
                                fn[-6:] == '.fasta']

    # Haal de promotorregionen van de genen gerelateerd aan de pathway op.
    print "Fetching promotor regions of %i genes related to pathway %s..." % (
        len(pathway_related_promotor_ids), met_pathway_code)
    get_specific_promotor_regions.main(pathway_related_promotor_ids,
                                       pathway_promotors_file,
                                       all_promotors_folder)

    # Scan alle motifs over de promotor regionen van de pathway.
    pathway_promotor_regions = list(
            SeqIO.parse(pathway_promotors_file, "fasta"))
    print "Scanning all %i motifs over all %i promotor regions related to" \
          " pathway %s..." % (
              len(all_human_motif_ids), len(pathway_promotor_regions),
              met_pathway_code)
    for pathway_promotor in pathway_promotor_regions:
        motif_scanner.scan_promotor(str(pathway_promotor.seq), pseudocount,
                                    pvalue, all_motifs_folder, all_motif_ids,
                                    pathway_promotor.id,
                                    pathway_scan_results_folder)

    analyze_scan_results.pandalize_csvs(all_motif_ids,
                                        pathway_scan_results_folder,
                                        "results_pathway_promotor_regions.csv")

    # Bepaal welke motif uniform voorkomt/hit in de set van promotorregionen
    #  gerelateerd aan de pathway.
    print "Calculating motif with highest coverage over all pathway" \
          " promotor regions..."
    pathway_motif_id = analyze_scan_results.get_best_motif(
            pathway_scan_results_folder + "merged_csv/",
            "results_pathway_promotor_regions.csv")

    for promotors_file in all_promotors_file_names:
        chromsosome_name = promotors_file.split("_")[-1][:-6]

        chromosome_promotors = list(
                SeqIO.parse(all_promotors_folder + promotors_file, "fasta"))
        print "Scanning best pathway-related motif over %i promotor regions " \
              "located on chromosome '%s' using %i CPU-cores..." % (
                  len(chromosome_promotors), chromsosome_name, CPU_count)

        motif_scanner.scan_chromosome(chromosome_promotors, pseudocount, pvalue,
                                      all_motifs_folder, pathway_motif_id,
                                      whole_genome_results_folder, CPU_count)

    # Verwerk de resultaten van de scan zodat er gerekend kan worden met de
    #  resultaten.
    analyze_scan_results.pandalize_csvs(pathway_motif_id,
                                        whole_genome_results_folder,
                                        "results_genome_promotor_regions.csv")

    print "Calculating top 100 high scoring promotor regions..."
    # Sorteer de score van de promotorregionen en bepaal de 100 hoogste.
    top_100_promotors = analyze_scan_results.get_top_100_scoring_genes(
            whole_genome_results_folder + "merged_csv/",
            "results_genome_promotor_regions.csv")

    # Haal de daadwerkelijke genen op en sla deze op in een multiple-fasta
    # bestand.
    print "Getting the genes of the top 100 high scoring promotor regions..."
    get_specific_genes.main(chromosomes_folder, chromosome_file_names,
                            top_100_promotors, top_100_genes_file)

    # Bepaal welke van de genen voldoen aan de eisen om deze daadwerkelijk
    # "coregulerende" genen te noemen.
    # Genen moeten minimaal 70 hits hebben.
    coregulated_genes = analyze_scan_results.get_coregulated_genes(
            whole_genome_results_folder + "merged_csv/",
            "results_genome_promotor_regions.csv", 70)

    # Split de top 100 in 50 en 20 voor de opvolgende biologische vragen.
    top_50_genes_file = get_top_x_genes.main(top_100_genes_file,
                                             top_100_genes_folder, 50)
    top_20_genes_file = get_top_x_genes.main(top_100_genes_file,
                                             top_100_genes_folder, 20)

    # Maak met ClustalW een multiple sequence alignment om de geconserveerde
    # regio vast te stellen.
    get_msa.main(top_50_genes_file, msa_file)
    # Bereken welk percentage van de ge-alignde genen geconserveerd is.
    percentage_conserved = get_percentage_conserved.main(msa_file)
    # Bepaal de naam van het langste gen van de 50 best scorende genen.
    longest_gene_name = get_name_longest_gene.main(top_50_genes_file)

    # Transleer de top 10 best scorende genen naar eiwitten. Extra genen worden
    # meegegeven indien 1 of meer van de 10
    # geen CDS heeft.
    translate_dna.main(top_100_genes_file, proteins_file, chromosomes_folder,
                       chromosome_file_names)

    # Bereken de phylogenetische boom van de 50 best scorende genen.
    get_paralogs_tree.main(top_50_genes_file, paralogs_file, phylo_tree_pic)

    # Bepaal van de top 20 genen wat de intron- en exon verhoudingen zijn en
    #  zet deze in een scatterplot.
    get_exon_sites.main(chromosomes_folder, chromosome_file_names,
                        top_20_genes_file, intron_exon_scatter_plot_pic)

    # Bepaal van de top 10 genen getransleert naar eiwitten de verhoudingen van
    # ladingen en zet deze in een frequentietabel.
    get_protein_charge_distribution.main(proteins_file,
                                         proteins_charge_distribution_pic)

    # Sla alle resultaten op in een Python Pickel object om later in te kunnen
    #  laden bij het genereren van de PDF.
    # Dit heeft als voordeel dat niet elke keer de pipeline opnieuw hoeft te
    #  draaien om een nieuwe render van de PDF te maken.
    resultaten_list = [
        # Auteurs
        auteurs,

        # studentnummer
        studentnummer,

        # gekozen metabolische route
        met_pathway_code,

        # De gekozen metabolische route als plaatje
        met_pathway_pic,

        # De transcriptiefactor die de gekozen metabolische route reguleert/
        # reguleren
        pathway_motif_id[0],

        # Het maximaal te vinden medegereguleerde genen buiten de
        # metabolische route
        coregulated_genes,

        # De naam van het langste gen
        longest_gene_name,

        # Het percentage van de geconserveerde regio van het langste gen
        percentage_conserved,

        # Een staafdiagram wordt gegenereerd.
        proteins_charge_distribution_pic,

        # Er wordt een fylogenetische boom gegenereerd
        phylo_tree_pic,

        # scatterplot
        intron_exon_scatter_plot_pic]

    pickle_file_object = open(pickle_file, "wb")
    pickle.dump(resultaten_list, pickle_file_object)
    pickle_file_object.close()

    # Genereer een PDF met de resultaten.
    pdf_generator.main(pickle_file)


if __name__ == "__main__":
    main(sys.argv)
