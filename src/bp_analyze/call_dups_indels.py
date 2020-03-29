from bg.tree import BGTree
from bg.multicolor import Multicolor
from bg.genome import BGGenome

from src.bp_analyze.tree_drawer import TreeDrawer
from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks
from src.bp_analyze.consistency_checker import TreeConsistencyChecker

import os

import matplotlib.pyplot as plt
import networkx as nx

sp_folder = 'data/E_coli/'
csv_file = sp_folder + 'total_stats.csv'
scale = 15000
blocks_folder = sp_folder + 'sibeliaz_out/fine/1000/'
grimm_file = blocks_folder + 'genomes_permutations.txt'
tree_file = sp_folder + 'tree/tree_midpoint_converted.nwk'
output_folder = blocks_folder + 'bp_results_dup_indels/'

if __name__ == "__main__":
    labels_dict = make_labels_dict(csv_file)
    tree_drawer = TreeDrawer(tree_file, scale=scale, labels_dict=labels_dict)
    consistency_checker = TreeConsistencyChecker(tree_file)

    genomes, blocks, block_genome_count = get_genomes_contain_blocks(grimm_file)
    print_species_stats(genomes, tree_drawer.get_all_leafs())

    os.makedirs(output_folder, exist_ok=True)
    for i, block in enumerate(blocks):
        print(i, 'of', len(blocks), ', current block:', block)

        genome_count = block_genome_count[block]
        labels = [f'{i}{"+" if i == tree_drawer.max_color - 1 else ""} copies' for i in range(max(genome_count.values()) + 1)]
        # print(labels)

        consistent = len(consistency_checker.check_consistency({genome: genome_count[genome] for genome in genomes})) == 0
        tree_drawer.draw(output_folder + f'block_{block}_{"" if consistent else "in"}consistent.pdf', colors_dict=genome_count, legend_labels=labels)
