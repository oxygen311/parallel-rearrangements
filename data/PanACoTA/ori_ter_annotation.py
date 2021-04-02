from iRep.gc_skew import gc_skew, parse_genomes
from datetime import datetime

import multiprocessing as mp

import pandas as pd
import numpy as np

import gzip
import json
import os.path
import time


window = 1000
slide = 10
plot_skew = False
workers = 32

gembase_file = 'Streptococcus_pyogenes/2-annotate_module/LSTINFO-LSTINFO-1314-filtered-1e-06_0.06.lst'
output_file = 'Streptococcus_pyogenes/ori_ters_test.csv'

gembase_df = pd.read_csv(gembase_file, sep='\t')
all_contigs = []


def run_subtask(sub_df, thread_id=0):
    ori_ters = []
    for i, (_, r) in enumerate(sub_df.iterrows()):
        print(f'{datetime.now()}: <thread {thread_id}>: {i} of {len(sub_df)}')

        with open(r.orig_name, 'r') as fasta_file:
            for name, length, seq in parse_genomes([fasta_file], False):
                # ori, ter, _1, _2 = gc_skew(name, length, seq, window, slide, plot_skew)
                ori, ter = 1, 2
                ori_ters.append([r.gembase_name, ori, ter])
                break # taking only first contig
    return ori_ters


if __name__ == '__main__':
    pool = mp.Pool()

    res = [pool.apply_async(run_subtask, [sub_df, thread_id])
           for thread_id, sub_df in enumerate(np.array_split(gembase_df, workers))]

    all_ori_ters = []
    for r in res:
        all_ori_ters.extend(r.get())

    df = pd.DataFrame(all_ori_ters, columns=['gembase_name', 'origin', 'terminus'])
    df.to_csv(output_file, sep=',', index=False)