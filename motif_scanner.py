#!/usr/bin/python
"""
Jeroen Merks
Multi-processing motif scanner. Op basis van MOODS.
J. Korhonen, P. Martinmaki, C. Pizzi, P. Rastas and E. Ukkonen. MOODS: fast
search for position weight matrix matches
in DNA sequences. Bioinformatics 25(23), pages 3181-3182. (2009)
C. Pizzi, P. Rastas and E. Ukkonen: Finding Significant Matches of Position Weight Matrices in Linear Time. IEEE/ACM
Transactions on Computational Biology and Bioinformatics. 8(1), pages 69 - 79. (2011)
"""
import Queue
import threading
from multiprocessing import Process

import MOODS.MOODS.parsers
import MOODS.MOODS.scan
import MOODS.MOODS.tools


def scan_promotor(seq, pseudocount, pvalue, all_motifs_folder, motif_ids,
                  promotor_id, results_folder):
    results_file = open(results_folder + "temp/" + promotor_id + ".csv", "w+")

    background = MOODS.MOODS.tools.bg_from_sequence_dna(seq, 1)

    matrices = [MOODS.MOODS.parsers.pfm_log_odds(all_motifs_folder + filename,
                                                 background, pseudocount) for
                filename in motif_ids]
    # reverse complements
    matrices = matrices + [
        MOODS.MOODS.parsers.pfm_log_odds_rc(all_motifs_folder + filename,
                                            background, pseudocount) for
        filename in motif_ids]
    thresholds = [MOODS.MOODS.tools.threshold_from_p(matrix, background, pvalue)
                  for matrix in matrices]
    # ---- end process all motifs ----

    # ---- scanning ----
    results = MOODS.MOODS.scan.scan_dna(seq, matrices, background, thresholds,
                                        7)
    # ---- end scanning ----

    # ---- process results ----
    # separate reverse complements and the non-reverse complements
    forward_r = results[:len(motif_ids)]
    reverse_r = results[len(motif_ids):]

    # mix the results together, use + and - to indicate strand
    results = [
        [(r.pos, r.score, '+') for r in forward_r[i]] + [(r.pos, r.score, '-')
                                                         for r in reverse_r[i]]
        for i in xrange(len(motif_ids))]

    results_file.write(promotor_id + ",")
    for result_nr, result in enumerate(results):
        if result_nr + 1 == len(results):
            results_file.write(str(len(result)))
        else:
            results_file.write(str(len(result)) + ",")

    results_file.close()
    # ---- end process results ----


def scan_chromosome(chromosome_promotors, pseudocount, pvalue,
                    all_motifs_folder, motif_ids, results_folder, cpu_count):
    # print "Queueing all promotor regions of chromosome '%s' for multi-processed scanning..." % chromosome_name
    queue = Queue.Queue()
    for promotor in chromosome_promotors:
        queue.put(promotor)

    def worker():
        while True:
            promotor = queue.get()
            # execute a task: call a shell program and wait until it completes
            process = Process(target=scan_promotor, args=(
                str(promotor.seq), pseudocount, pvalue, all_motifs_folder,
                motif_ids, str(promotor.id), results_folder))
            process.start()
            process.join()  # this blocks until the process terminates
            queue.task_done()

    for cpu in range(cpu_count):
        thread = threading.Thread(target=worker)
        thread.daemon = True
        thread.start()
    queue.join()  # block until all tasks are done
