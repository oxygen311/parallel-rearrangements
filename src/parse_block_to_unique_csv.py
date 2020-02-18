from io import StringIO
from collections import Counter

import numpy as np
import pandas as pd

import csv

dir = 'sibeliaz/Streptococcus_suis/fine/5000/'
file = 'blocks_coords.txt'
out_file = 'blocks_coords_unique_gene.csv'
sep = '-' * 80

if __name__ == "__main__":
    lines = open(dir + file).read().split('\n')[:-1]
    ls = np.array(lines)
    bs = np.split(ls, np.where(ls == sep)[0])

    # names of chromosomes
    df_names = pd.read_csv(StringIO('\n'.join(bs[0])), sep='\t')
    seq_names = {}
    for index, row in df_names.iterrows():
        seq_names[row['Seq_id']] = row['Description']

    # blocks data
    len_seqs = len(seq_names.keys())
    blocks_count = 0

    ans = [['Block', 'Seq', 'Start', 'End', 'Strand']]
    for b in bs[1:-1]:
        block = b[1].split('#')[1]
        df_block = pd.read_csv(StringIO('\n'.join(b[2:])), sep='\t')

        seq_ids = df_block['Seq_id'].values
        cnt = Counter(seq_ids)

        if not all(c == 1 for c in cnt.values()): continue
        if len(set(seq_ids)) != len_seqs: continue

        blocks_count += 1
        for index, row in df_block.iterrows():
            seq_name, start, end, strand = seq_names[row['Seq_id']], row['Start'], row['End'], row['Strand']
            ans.append([block, seq_name, start, end, strand])

    with open(dir + out_file, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(ans)

    print(f'Wrote {blocks_count} blocks of {len_seqs} sequences in {len(ans) - 1} lines')