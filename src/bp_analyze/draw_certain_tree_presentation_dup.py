from bg.tree import BGTree
from bg.multicolor import Multicolor
from bg.genome import BGGenome

from src.bp_analyze.tree_holder import TreeHolder
from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks
from src.bp_analyze.shigella_common import get_class_colors, count_shigella_differs_from_value

from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import MDS
from collections import defaultdict

import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sp_folder = 'data/E_coli/'
csv_file = sp_folder + 'total_stats.csv'
scale = 15000
blocks_folder = sp_folder + 'sibeliaz_out/fine/1000/'
grimm_file = blocks_folder + 'genomes_permutations.txt'
tree_file = sp_folder + 'tree/tree_midpoint_converted.nwk'
output_folder = blocks_folder + 'supplementary/copy_number_variation/old_version/'


if __name__ == "__main__":
    # labels_dict = make_labels_dict(csv_file)

    # consistency_checker = TreeConsistencyChecker(tree_file)

    genomes, blocks, block_genome_count = get_genomes_contain_blocks(grimm_file)
    # print_species_stats(genomes, tree_drawer.get_all_leafs())

    for block, count in block_genome_count.items():
        if any(c > 1 for c in count.values()):
            print(block)

    # interesting_blocks = ['1368', '910', '911', '912', '919', '920', '270', '271', '908', '1362', '1364', '1360']
    # for block in interesting_blocks:
    #
    #     count = block_genome_count[block]
    #
    #     tree_drawer = TreeHolder(tree_file, scale=scale, labels_dict=labels_dict)
    #     labels = [f'{i}{"+" if i == tree_drawer.max_color - 1 else ""} copies'
    #               for i in range(max(count.values()) + 1)]
    #     tree_drawer.draw(output_folder + f'block_{block}.pdf',
    #                      colors_dict=count, show_branch_support=True,
    #                      show_scale=True, legend_scale=1, draw_pres=True)

