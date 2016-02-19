"""
Christian van der Niet
Programma maakt op basis van een gegeven multiple-fasta DNA bestand een
multiple sequence alignment.

"""
import os


def main(top_50_in_file, alignment_out_file):
    # de wrapper in biopython werkt niet op elk systeem.
    # hierdoor word clustalw direct ->
    # aangeroepen met een command line comando
    # functie heeft input file nodig in fasta formaat
    # deze word zonder extencie doorgegeven
    print "Generating multiple sequence alignment..."
    os.system("clustalw -INFILE=%s -OUTFILE=%s -QUIET" % (
    top_50_in_file, alignment_out_file))
