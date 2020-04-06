from bg.tree import BGTree
from bg.multicolor import Multicolor
from bg.genome import BGGenome

from src.bp_analyze.tree_drawer import TreeDrawer
from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks
from src.bp_analyze.consistency_checker import TreeConsistencyChecker
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
output_folder = blocks_folder + 'bp_results_copy_number_variation/'
threshold = 10


def construct_dist_matrix(counts):
    n = len(counts)
    arr = np.zeros((n, n))
    for i1, count1 in enumerate(counts):
        for i2, count2 in enumerate(counts):
            if i1 >= i2: continue
            arr[i1, i2] = arr[i2, i1] = \
                sum(count1[genome] != count2[genome] for genome in count1.keys() | count2.keys())

    return arr


if __name__ == "__main__":
    labels_dict = make_labels_dict(csv_file)
    tree_drawer = TreeDrawer(tree_file, scale=scale, labels_dict=labels_dict)
    consistency_checker = TreeConsistencyChecker(tree_file)

    genomes, blocks, block_genome_count = get_genomes_contain_blocks(grimm_file)
    print_species_stats(genomes, tree_drawer.get_all_leafs())

    # os.makedirs(output_folder, exist_ok=True)
    possible_counts = [block_genome_count[block] for block in blocks]
    arr = construct_dist_matrix(possible_counts)
    print("<Matrix done>")
    cls = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average', distance_threshold=threshold).fit_predict(arr)
    print('<Clustering done>')
    print('Clusters:', np.unique(cls).shape[0])

    sns.set(style="whitegrid", font="serif", font_scale=0.3)
    for cl in np.unique(cls):
        inds = np.where(cls == cl)[0]
        print(inds.shape)
        cur_blocks = [blocks[i] for i in inds]
        cur_counts = [possible_counts[i] for i in inds]

        inconsistent = False
        max_dups, max_inds = 0, 0
        count_blocks = defaultdict(list)
        for block, count in zip(cur_blocks, cur_counts):
            class_colors = get_class_colors(labels_dict, {genome: count[genome] for genome in genomes})
            e_coli_mean = np.mean(class_colors['Escherichia coli'])
            shigella_differs_dup, shigella_differs_indel = count_shigella_differs_from_value(class_colors, e_coli_mean)

            inconsistent |= len(consistency_checker.check_consistency({genome: count[genome] for genome in genomes})) != 0
            max_dups = max(max_dups, shigella_differs_dup)
            max_inds = max(max_inds, shigella_differs_indel)
            count_blocks[frozenset(count.items())].append(block)

        type = 'none'
        if max_dups > 0 and max_inds > 0: type = 'both'
        elif max_dups > 0: type = 'duplication'
        elif max_inds > 0: type = 'indel'

        curr_output_folder = output_folder + type + '/' \
                             + f'{"" if not inconsistent else "in"}consistent__{max_dups + max_inds}_shigells_differ/' \
                             + f'cluster_{cl}_size_{len(cur_blocks)}/'
        os.makedirs(curr_output_folder, exist_ok=True)

        for i, (unique_count, unique_count_blocks) in enumerate(count_blocks.items()):
            count = defaultdict(int, unique_count)
            labels = [f'{i}{"+" if i == tree_drawer.max_color - 1 else ""} copies'
                      for i in range(max(count.values()) + 1)]
            tree_drawer.draw(curr_output_folder + f'tree_{i}_representing_{len(unique_count_blocks)}_blocks.pdf',
                             colors_dict=count, legend_labels=labels)

            with open(curr_output_folder + f'tree_{i}_representing_{len(unique_count_blocks)}_blocks.txt', 'w') as f:
                print('\n'.join(map(str, unique_count_blocks)), file=f)
