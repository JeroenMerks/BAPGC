#!/usr/bin/python
"""
Jeroen Merks
Download en verwerk een gegeven lijst van motifnamen van de meest recente versie
 van de JASPAR CORE database.

"""
import urllib


def main(motif_ids, motifs_folder):
    print 'Downloading and processing all %i motifs of transcription factor ' \
          'binding sites known to Homo sapien from current JASPAR_CORE' % len(
        motif_ids)

    for motif_id in motif_ids:
        motif_file_name = motif_id + ".pfm"
        pfm_file = open(motifs_folder + motif_file_name, "w+")
        fetched_text_file = urllib.urlopen(
                "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/" + motif_file_name).readlines()
        for lin_nr, line in enumerate(fetched_text_file[1:]):
            if lin_nr == 3:
                pfm_file.write(line[5:-3])
            else:
                pfm_file.write(line[5:-3] + "\n")

        pfm_file.close()
