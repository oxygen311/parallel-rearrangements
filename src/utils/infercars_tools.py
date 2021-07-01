import numpy as np
import re
import pandas as pd

from io import StringIO


p = "([A-Za-z0-9_\(\)\/\s]+)\.([A-Za-z0-9_\.]+):(\d+)-(\d+) ([+|-]).*"
pattern = re.compile(p)
columns = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]


def find_indices(lst, condition):
    return [i for i, elem in enumerate(lst) if condition(elem)]

def parse_to_df(file_name):
    with open(file_name) as f:
        lines = f.readlines()

    last_line = len(lines) - 1
    while lines[last_line] == '\n': last_line -= 1

    n_at_end = len(lines) - 1 - last_line
    for _ in range(1 - n_at_end): lines.append('\n')

    bs = np.split(lines, find_indices(lines, lambda x: x[0] == ">"))
    temp = []

    # b_inds = [int(lines[i][1:]) for i in inds]

    for i, b in enumerate(bs):
        if len(b) == 0: continue
        b_i = int(b[0][1:])

        for oc in b[1:-1]:
            m = pattern.match(oc)
            temp.append([b_i, m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)])

    return pd.DataFrame(temp, columns=columns)

def export_df(df, file_name):
    with open(file_name, 'w') as f:
        for block, block_df in df.groupby('block'):
            print(f'>{block}', file=f)
            for i, row in block_df.iterrows():
                print(f'{row["species"]}.{row["chr"]}:{row["chr_beg"]}-{row["chr_end"]} {row["orientation"]}', file=f)
            print(file=f)

def filter_unique_gene(df):
    allowed_blocks = set()
    all_sp = len(df['species'].unique())
    for block, df_block in df.groupby('block'):
        if len(df_block) == len(df_block['species'].unique()) == all_sp:
            allowed_blocks.add(block)

    return df[df.apply(lambda x: x['block'] in allowed_blocks, axis=1)]

def dist_between_blocks(df):
    ds = []
    for sp, df_sp in df.groupby('species'):
        df_sp = df_sp.sort_values(by=['chr_beg'])
        ds += (start_ - end_ for start_, end_ in zip(df_sp['chr_beg'][1:], df_sp['chr_end']))
    return ds

def blocks_length_dist(df):
    return [end_ - start_ for start_, end_ in zip(df['chr_beg'], df['chr_end'])]

def number_of_genomes_dist(df):
    return [len(df_block.species.unique()) for _, df_block in df.groupby('block')]


def ragout_to_infercars(dir, in_file, out_file):
    sep = '-' * 80
    split_by_underscore = False

    lines = open(dir + in_file).read().split('\n')[:-1]
    ls = np.array(lines)
    bs = np.split(ls, np.where(ls == sep)[0])

    # names of chromosomes
    df_names = pd.read_csv(StringIO('\n'.join(bs[0])), sep='\t')
    chr_names = {}
    for index, row in df_names.iterrows():
        chr_names[row['Seq_id']] = row['Description']

    # blocks data
    with open(dir + out_file, 'w') as f:
        for b in bs[1:-1]:
            block = b[1].split('#')[1]
            df_block = pd.read_csv(StringIO('\n'.join(b[2:])), sep='\t')

            print(f'>{block}', file=f)
            for index, row in df_block.iterrows():
                chr_name, start, end, strand = chr_names[row['Seq_id']], row['Start'], row['End'], row['Strand']
                if split_by_underscore:
                    splitted = chr_name.split('_', 1)
                    chr = splitted[0] + '.' + splitted[1].split('.')[0]
                else:
                    chr = chr_name
                if not 'alt' in chr:
                    print(chr + ':' + (f'{start}-{end}' if strand == '+' else f'{end}-{start}') + ' ' + strand, file=f)
            print(file=f)