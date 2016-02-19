#!/usr/bin/python
# -*- coding: utf-8 -*-

# Melissa van Wieringen
# s1079422
# Python PDF generator
# Last changes: 2 februari 2016


# Imports
import pickle
import re
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, \
    Image


# Functie die een Pickle Object uitpakt en returnt.
def unpack_pickle(pickle_object_location):
    pickle_file_object = open(pickle_object_location, "rb")
    results_list = pickle.load(pickle_file_object)
    return results_list


# Class generator om pdf bestanden mee te creeeren.
class generator:
    # Init.
    def __init__(self):

        # Maak een nieuw document
        self.doc = SimpleDocTemplate("BAPGC.pdf", pagesize=letter,
                                     rightMargin=72, leftMargin=72,
                                     topMargin=72, bottomMargin=18)
        # Maak een nieuwe story list.
        self.story = []

        # Maak een nieuw StyleSheet om styles aan toe te voegen
        self.styles = getSampleStyleSheet()

        self.styles.add(ParagraphStyle(name='Header1', fontSize=18, leading=24,
                                       spaceAfter=10))

        self.styles.add(ParagraphStyle(name='Header2', fontSize=16, leading=24,
                                       spaceAfter=10))

        self.styles.add(ParagraphStyle(name='Header3', fontSize=14, leading=24,
                                       spaceAfter=10))

    # add_pagenumber. Ontvangt het canvas en het document en zet een pagina
    # nummer onderin de pagina
    def add_pagenumber(self, canvas, doc):

        # Krijgt het paginanummer terug
        pagina_nummer = canvas.getPageNumber()

        text = "Pagina %s" % pagina_nummer

        # Tekent het nummer rechtsonderin de pagina
        canvas.drawRightString(200 * mm, 20 * mm, text)

    # generate_page. Creeert een pagina. Ontvangt een chapter list.
    # Zoekt naar de volgende keywords:
    # __SPACER_y_xxx    waar y de width is en xxx de height is.
    #                   (y is altijd 1 getal!)
    # __PICTURE_x_y__zz waar x de heigth is, y de width
    #                   en zz de url is (in de huidige dir!)
    #                   Opmerking: plaats altijd twee underscores (_) tussen
    #                   de width en de url!!
    # __HEADER_y_xx     waar y 1, 2 of 3 is
    #                   1 = grote header; 2 = middel header; 3 = kleine header)
    #                   en xx de string is
    def generate_page(self, chapter):

        # Controleer of chapter wel inhoud heeft en of het wel een lijst is
        if len(chapter) <= 0 or type(chapter) != list:
            print "FOUTMELDING: Je geeft een lege pagina aan generate_page."
            return False

        # Append de inhoud lijst aan de story, maar sla de eerste index
        # (de parameters) en de tweede index (titel) over.
        for paragraph in chapter:

            if "__SPACER" in paragraph:
                # Geef een spacer mee
                self.story.append(
                        Spacer(float(paragraph[9]), int(paragraph[11:])))

                # Zorg dat deze parameter niet geprint wordt
                continue

            elif "__HEADER" in paragraph:

                if int(paragraph[9]) in [1, 2, 3]:
                    temp_var = "Header" + paragraph[9]
                    # Geef deze regel een groot letterype
                    self.story.append(
                            Paragraph(paragraph[11:], self.styles[temp_var]))

                else:
                    print "FOUTMELDING: Je hebt geen 1, 2 of 3 op positie 9 " \
                          "van " + paragraph + " staan!"
                    return False

                # Zorg dat deze parameter niet geprint wordt
                continue

            elif "__PICTURE" in paragraph:

                # Alles behalve __PICTURE_
                x = paragraph[10:]

                # Alles gesplit op de _
                splitted = x.split("_")

                width = float(splitted[0])
                height = float(splitted[1])

                # Het splitten op __ om de naam van het plaatje te achterhalen
                name = re.compile("__").split(paragraph)
                name = name[-1]

                # Maakt een nieuw image object aan met de meegegeven url
                # en width en height
                im = Image(name, width * mm, height * mm)
                self.story.append(im)

                # Zorg dat deze parameter niet geprint wordt
                continue

            self.story.append(Paragraph(paragraph, self.styles["Normal"]))

        # Na elke pagina (chapter) komt een pageBreak om een nieuwe pagina
        # te beginnen
        self.story.append(PageBreak())

    # generate_pdf. Creeert de pdf met de story.
    def generate_pdf(self):

        # Bouw de pdf met de story.
        self.doc.build(self.story, onLaterPages=self.add_pagenumber)


# main. Bouwt de pdf met behulp van de generator class.
def main(pickle_file_location):
    results = unpack_pickle(pickle_file_location)

    # results = [
    #     'Melissa van Wieringen, Jeroen Merks, Koen van Diemen, Rick de Graaf en Christian van der Niet',
    #     's1072614', 'hsa04915', 'data/results/met_pathway_picture/hsa04915.png',
    #     'MA0481.1.pfm',
    #     ['LOC105376066', 'DAPK3', 'SMARCA4', 'VMAC', 'DNASE2', 'MINK1',
    #      'C17orf53', 'ARMC6', 'LOC101928572', 'COL26A1', 'DAND5', 'ZSWIM4',
    #      'FTL', 'PRR12', 'MIR3188', 'FAHD1', 'KIF1C', 'ARRDC1', 'RNF207',
    #      'PUS1', 'VPS37D', 'VARS', 'SCAMP4', 'KRI1', 'GIPR', 'TALDO1', 'CARNS1',
    #      'RPL23AP5', 'LOC105376712', 'DCAF15', 'SLC12A4', 'CCDC151', 'PDLIM7',
    #      'ADAT3', 'RASA4CP', 'MIR4516', 'CENPM', 'LOC100419925', 'ZNF414',
    #      'BTBD2', 'TJP1P', 'BCKDK', 'RPS15P9', 'LOC101928543', 'PKN1', 'ZNF487',
    #      'LTC4S', 'ATAD3B', 'MAP2K2', 'SELV', 'ZNF101', 'LOC100507373', 'CDC34',
    #      'RPL32P34', 'DGCR8', 'MIR5090', 'MIR1281', 'LOC100129352', 'KLF1',
    #      'CORO1A', 'FBN3', 'EP300', 'MIR1199', 'MIR6798', 'C19orf43', 'PPP1R26',
    #      'LOC100419924', 'PLP2', 'SPG7', 'APOBEC3D', 'SPIRE2', 'APBA3',
    #      'CATSPERD', 'LOC105370690', 'TSSK6', 'C16orf90', 'LOC105372266',
    #      'APC2', 'NDUFA13', 'ZDHHC12', 'ATP1A3', 'GALR3', 'MAPK8IP2', 'CARM1',
    #      'GMIP'], 'SMARCA4', 0.0,
    #     'data/results/proteins/protein_charge_distribution.png',
    #     'data/results/phylo_tree/phylogenetic_tree_paralogs.png',
    #     'data/results/intron_exon_scatter_plot/intron_exon_scatter_plot.png']

    ########################################
    # HIER KOMEN DE PAGINA'S VOOR IN DE PDF#
    ########################################
    front_page = [
        "__HEADER_1_Onderzoek coregulatie genen buiten metabolische route:",
        "__HEADER_2_" + results[2], "__PICTURE_150_150__" + results[3],
        "KEGG representatie van de metabole route " + results[2],
        "__SPACER_1_100", "Auteurs: " + results[0],
        "Studentnummer: " + results[1], "Datum: Februari 2016"]
    second_page = ["__HEADER_1_Inleiding",
                   "Dit is een automatisch gegenereerd rapport van een "
                   "analyse  op coregulatie van genen buiten de volgende "
                   "metabolische route: " + results[2] + ".",
                   "De conclusie en discussie van de resultaten zijn mogelijk "
                   "later met de hand toegevoegd.", "__SPACER_1_10",
                   "Het doel van de pipeline is het vinden van genen buiten de"
                   " pathway die door eenzelfde transcriptiefactor worden "
                   "gereguleerd.",
                   "Aan de hand van 5 onderzoeksvragen worden deze "
                   "medegereguleerde genen verder onderzocht.",
                   "Achtergrondinformatie over de gekozen metabolische route", ]
    third_page = ["__HEADER_1_Materialen en methoden",
                  "114 Motifs van transcriptiefactoren gerelateerd aan de "
                  "mens, zoals bekend bij de laatste JASPAR CORE datebase ("
                  "minus motifs MA0528.1 en MA0050.2 vanwege overrepresentatie "
                  "in eerder gedane tests), ",
                  "zijn met een p-value van 0.0001 gescant over álle "
                  "promotorregionen (6000 basenparen vóór de start van het "
                  "gen en 500 ná) van alle genen die volgens KEGG gerelateert "
                  "zijn aan de opgegeven metabole route.",
                  "Er is ervoor gekozen om de beste transcriptiefactor te "
                  "selecteren om over de promotorregionen van de rest van het "
                  "genoom te scannen. Dit om het onderzoek duidelijk af te "
                  "kaderen en er met minder ruis in de resultaten betere "
                  "conclusies te trekken zijn op mogelijke corelaties tussen "
                  "de pathway en gevonden genen waar motifs van de "
                  "transcriptiefactoren sterk op hitten.",
                  "Het criterium voor het bepalen van de beste coregulerende "
                  "transcriptiefactor is een zo groot mogelijke overlap van "
                  "hits over alle promotorregionen van de pathway.",
                  "Indien meerdere transcriptiefactoren een coverage van 100%"
                  " bleken te hebben, werd de transcriptiefactor met "
                  "cumulatief het hoogste aantal hits genomen.",
                  "__SPACER_1_10",
                  "De resultaten van de scan op het hele genoom zijn "
                  "vervolgens stringent geselecteerd, zodat een handje vol "
                  "genen over zouden blijven. Hier is voor gekozen zodat "
                  "vervolgonderzoek kan worden gedaan op basis van resultaten "
                  "die sterk uit de analyse naar voren zijn gekomen en het "
                  "waarschijnlijk waard zijn om met de hand verder te "
                  "onderzoeken.", "__SPACER_1_10",
                  "Bij een p-value van 0.0001 en een minimum aantal hits "
                  "van 100 van de transcriptiefactor in kwestie, kwam dit "
                  "neer op " + str(len(results[5])) + " promotors."]
    fourth_page = ["__HEADER_1_Resultaten",
                   "De volgende motif van transcriptiefactor " + results[
                       4] + " had de hoogste coverage over de "
                            "promotorregionen van de pathway.",
                   "De gensymbolen van de volgende genen hadden minimaal 70 "
                   "hits van deze transcriptiefactor op de promotorregio: " + str(
                           results[5]), "__SPACER_1_10",
                   "Van de 50 beste medegereguleerde genen is een local "
                   "multiple sequence alignment gemaakt met behulp van "
                   "ClustalW 2.1, gebruikmakende van de standaard parameters.",
                   "Alle genen korter dan gen: " + results[
                       6] + " vormden een geconserveerde regio van" + str(
                           results[7]) + "% ten opzichte van dat gen.",
                   "__SPACER_1_10",
                   "De 10 beste medegereguleerde genen zijn getransleert naar "
                   "hun respectievelijke eiwit met de SeqIO.translate("
                   "cds=True) functie van BioPython 1.66.",
                   "Niet van elk gen dat in de pipeline was onderzocht was de"
                   " CDS (CoDing Sequence) bekend. Indien dit het geval was "
                   "dan werd het eerstvolgende "
                   "(minder hoor scorende) gen uit de gesorteerde lijst van "
                   "hoog scorende genen geselecteerd.", "__SPACER_1_10",
                   "Eiwitten bestaan uit aminozuren die elk of hydrofiel, "
                   "hydrofoob of neutraal zijn.",
                   "Van de 10 eiwitten is een staafdiagram, te zien in Fig. "
                   "1, gemaakt waar de verhoudingen van de hydrofobiciteit "
                   "van de aminozuren uit af te lezen is.",
                   "__PICTURE_120_80__" + results[8],
                   "Fig. 1, staafdiagram van de verhoudingen van de "
                   "hydrofobiciteit van de aminozuren van de 10 eiwitten.",
                   "__SPACER_1_10",
                   "Voor elk van de 10 beste medegereguleerde genen zijn de 4 "
                   "meest verwante paraloge genen bepaald.",
                   "Deze paraloge genen zijn vervolgens in een fylogenetische "
                   "boom, te zien in Fig. 2, uiteengezet.",
                   "__PICTURE_80_80__" + results[9],
                   "Fig. 2, de phylogenetische boom van de 4 meest verwante "
                   "paraloge genen per gen van de 10 beste genen, gemaakt met "
                   "ClustalW.", "__SPACER_1_10",
                   "Van de 20 beste medegereguleerde genen zijn de cumulatieve"
                   " intron- en exonlengte in een scatterplot, te zien in "
                   "Fig. 3, uiteen gezet.",
                   "Ook van deze genen was het niet altijd bekend waar de"
                   " exonen (en dus ook intronen) zich bevonden. Deze zijn "
                   "uit de scatterplot weggelaten.", "__SPACER_1_10",
                   "__PICTURE_150_120__" + results[10],
                   "Fig. 3, een scatterplot van de de cumulatieve intron- en "
                   "exonlengtes van de 20 beste genen."]
    fifth_page = ["__HEADER_1_Conclusie & Discussie"]
    sixth_page = ["__HEADER_1_Referenties"]

    print "Generating PDF of results..."

    # Maak een nieuwe generator.
    pdf_maker = generator()

    # Genereer al de pagina's, een voor een.
    # Op het einde wordt de pdf gecreeerd met al de pagina's.
    try:
        pdf_maker.generate_page(front_page)
        pdf_maker.generate_page(second_page)
        pdf_maker.generate_page(third_page)
        pdf_maker.generate_page(fourth_page)
        pdf_maker.generate_page(fifth_page)
        pdf_maker.generate_page(sixth_page)

        pdf_maker.generate_pdf()
    except IOError:
        print "FOUTMELDING: Je hebt ergens geen juiste bestandsnaam meegegeven."
    except AttributeError:
        print "FOUTMELDING: Attribuut error."
    except ValueError:
        print "FOUTMELDING PDF generator is gestopt door een fout."

    print "PDF generator is done working. Thank you for your cooperation and " \
          "come again!"

# main()
