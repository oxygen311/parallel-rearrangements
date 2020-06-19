import csv

import pandas as pd
import numpy as np

from collections import defaultdict, Counter

from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks, get_block_seqs

block_size = 250
folder = f'data/EPEC/blocks/{block_size}/'

grimm_file = 'genomes_permutations.txt'
csv_file = 'data/EPEC/aEPECmeta.csv'

# row_of_interest = 'inflammatory_potential_in_vitro'
row_of_interest = 'disease_healthy'
not_interesting = {'unknown'}

def get_number_to_type():
    number_to_type = {}
    for i, row in df.iterrows():
        number, type = row.EPEC_number, row[row_of_interest]
        if type in not_interesting: continue
        number_to_type[number] = type
    return number_to_type

def get_groups(count):
    d = defaultdict(list)
    for number, copies in count.items():
        d[copies].append(number)
    return list(d.values())

def information(group, number_to_type):
    types = [number_to_type[n] for n in group]
    cnt = Counter(types)
    print(cnt)

    ps = np.array([v / len(group) for v in cnt.values()])
    return - np.sum(ps * np.log2(ps))


if __name__ == "__main__":
    genomes, blocks, block_genome_count = get_genomes_contain_blocks(folder + grimm_file, skip_unique=True,
                                                                     split_by='.')
    df = pd.read_csv(csv_file)

    number_to_type = get_number_to_type()
    rows = [['block', 'gain']]

    ns = number_to_type.keys()
    all_information = information(ns, number_to_type)

    for block in blocks:
        count = {n: block_genome_count[block]['EPEC' + str(n)] for n in ns}
        groups = get_groups(count)

        information_gain = all_information

        for group in groups:
            group_i = information(group, number_to_type)
            information_gain -= len(group) / len(ns) * group_i

        print(block, information_gain)
        rows.append([block, information_gain])

    with open(folder + f'information_gain_{row_of_interest}.csv', 'w') as f:
        wtr = csv.writer(f)
        wtr.writerows(rows)

    print('Max gain:', max(row[1] for row in rows[1:]))
