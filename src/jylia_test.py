from bg.grimm import GRIMMReader
from bg.tree import BGTree
from bg.genome import BGGenome
from bg.vertices import BGVertex

from src.bp_analyze.common import get_genomes_contain_blocks, print_species_stats, make_labels_dict
from src.bp_analyze.common_bg import draw_bp_with_weights, get_colors_by_edge, save_pygraphviz
from src.bp_analyze.consistency_checker import TreeConsistencyChecker
from src.bp_analyze.shigella_common import get_class_colors, white_proportion, count_shigella_differs_from_white

from pprint import pprint

import pandas as pd

import os

sp_folder = 'data/Staphylococcus_aureus/'
csv_file = sp_folder + 'assemblies_chrs.csv'
tree_file = sp_folder + '/tree/RAxML_bipartitionsBranchLabels.attepmt_converted'

if __name__ == "__main__":
    df = pd.read_csv(csv_file)
    print(len(df.chromosome.unique()))

    with open(tree_file) as f:
        l = f.readline()

    for i, row in df.iterrows():
        chr, assembl = row['chromosome'], row['assembly']
        l = l.replace(assembl, chr)

    with open(tree_file + '_renamed', 'w') as f:
        print(l, file=f)