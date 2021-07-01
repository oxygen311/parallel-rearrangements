from bg.tree import BGTree
from bg.multicolor import Multicolor
from bg.genome import BGGenome

from src.bp_analyze.tree_holder import TreeHolder
from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks, get_block_seqs
from src.bp_analyze.shigella_common import get_class_colors, count_shigella_differs_from_value

from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import MDS
from collections import defaultdict
from itertools import chain
from pprint import pprint

import os
import csv

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sp_folder = 'data/Staphylococcus_aureus/'
csv_file = sp_folder + 'assemblies_chrs.csv'
scale = 72000
blocks_folder = sp_folder + 'sibeliaz/fine/1000/'
grimm_file = blocks_folder + 'genomes_permutations.txt'
tree_file = sp_folder + 'tree/RAxML_bipartitionsBranchLabels.attepmt_converted_renamed_midpoint'
output_folder = blocks_folder + 'copy_number_variation/'
csv_out_file = output_folder + 'stats.csv'
threshold = 0.125
j = 0.75
b = 0.25

def construct_dist_matrix_similarity(counts):
    n = len(counts)
    arr = np.zeros((n, n))
    for i1, count1 in enumerate(counts):
        for i2, count2 in enumerate(counts):
            if i1 >= i2: continue
            arr[i1, i2] = arr[i2, i1] = \
                sum(count1[genome] != count2[genome] for genome in count1.keys() | count2.keys())

    return arr


def construct_dist_matrix_proximity(blocks, seqs, percentile=25):
    def block_to_positions_index(seq):
        positions = defaultdict(list)
        for i, block in enumerate(seq):
            positions[block].append(i)
        return positions

    def get_distance(indexes, seqs, block1, block2):
        ds = []
        for ind, seq in zip(indexes, seqs):
            n = len(seq)
            seq_ds = [(abs(pos1 - pos2) + n) % n for pos1 in ind[block1] for pos2 in ind[block2]]
            if len(seq_ds) > 0:
                ds.append(min(seq_ds))
        if len(ds) == 0:
            return max_len_seq
        else:
            return np.percentile(ds, percentile)

    indexes = [block_to_positions_index(seq) for seq in seqs]
    max_len_seq = max(map(len, seqs))

    n = len(blocks)
    arr = np.zeros((n, n))
    for i1, block1 in enumerate(blocks):
        for i2, block2 in enumerate(blocks):
            if i1 >= i2: continue
            arr[i1, i2] = arr[i2, i1] = get_distance(indexes, seqs, block1, block2)

    return arr


if __name__ == "__main__":
    labels_dict = make_labels_dict(csv_file)
    tree_holder = TreeHolder(tree_file, scale=scale, labels_dict=labels_dict)

    genomes, blocks, block_genome_count = get_genomes_contain_blocks(grimm_file, skip_unique=True)
    print_species_stats(genomes, tree_holder.get_all_leafs())

    print(len(blocks))
    blocks = [block for block in blocks if any(v > 1 for v in block_genome_count[block].values())]
    print(len(blocks))

    # os.makedirs(output_folder, exist_ok=True)
    possible_counts = [block_genome_count[block] for block in blocks]
    J = construct_dist_matrix_similarity(possible_counts)
    # J = np.sqrt(J)
    J /= J.max()

    seqs = get_block_seqs(grimm_file)
    B = construct_dist_matrix_proximity(blocks, seqs)
    # B = np.sqrt(B)
    B /= B.max()

    D = J * j + B * b
    B /= D.max()

    print("<Matrix done>")
    cls = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average',
                                  distance_threshold=threshold).fit_predict(D)
    print('<Clustering done>')
    print('Blocks:', len(blocks))
    print('Clusters:', np.unique(cls).shape[0])

    ans = [['cluster', 'tree', 'parallel_rear_score', 'parallel_rear_unique_innovation_count',
            'parallel_rear_all_innovations_count']]

    sns.set(style="whitegrid", font="serif", font_scale=0.3)
    for cl in np.unique(cls):
        inds = np.where(cls == cl)[0]
        print(inds.shape)
        cur_blocks = [blocks[i] for i in inds]
        cur_counts = [possible_counts[i] for i in inds]

        count_blocks = defaultdict(list)
        for block, count in zip(cur_blocks, cur_counts):
            count_blocks[frozenset(count.items())].append(block)

        curr_output_folder = output_folder + f'cluster_{cl}_size_{len(cur_blocks)}/'
        os.makedirs(curr_output_folder, exist_ok=True)

        for i, (unique_count, unique_count_blocks) in enumerate(count_blocks.items()):
            count = defaultdict(int, unique_count)

            # consistency
            tree_holder.count_innovations_fitch({genome: count[genome] for genome in genomes})

            score_rear, count_rear, count_all_rear = tree_holder.count_parallel_rearrangements(skip_grey=False)

            ans.append([cl, i, score_rear, count_rear, count_all_rear])

            labels = [f'{i}{"+" if i == tree_holder.max_color - 1 else ""} copies'
                      for i in range(max(count.values()) + 1)]
            tree_holder.draw(curr_output_folder + f'tree_{i}_representing_{len(unique_count_blocks)}_blocks.pdf',
                             legend_labels=labels, show_branch_support=False, legend_scale=0.42)

            with open(curr_output_folder + f'tree_{i}_representing_{len(unique_count_blocks)}_blocks.txt', 'w') as f:
                print('\n'.join(map(str, unique_count_blocks)), file=f)

    with open(csv_out_file, 'w') as f:
        wtr = csv.writer(f)
        wtr.writerows(ans)