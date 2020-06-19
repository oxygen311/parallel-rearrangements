from bg.grimm import GRIMMReader
from bg.tree import BGTree
from bg.genome import BGGenome
from bg.vertices import BGVertex

from src.bp_analyze.tree_holder import TreeHolder
from src.bp_analyze.common import get_genomes_contain_blocks, print_species_stats, make_labels_dict
from src.bp_analyze.common_bg import draw_bp_with_weights, get_colors_by_edge, save_pygraphviz

import numpy as np

import os
import csv

sp_folder = 'data/Staphylococcus_aureus/'
csv_file = sp_folder + 'assemblies_chrs.csv'
scale = 72000
blocks_folder = sp_folder + 'sibeliaz/fine/1000/'
grimm_file = blocks_folder + 'genomes_permutations_unique.txt'
tree_file = sp_folder + 'tree/RAxML_bipartitionsBranchLabels.attepmt_converted_renamed_midpoint'
output_folder = blocks_folder + 'bp_unique_block_output/'
csv_out_file = output_folder + 'stats.csv'

def white_proportion(colors):
    return np.mean(list(map(lambda c: c == 0, colors)))

def get_genome_colors_by_edge(bg, edge, genomes):
    def get_neighbour_with_genome(v, genome):
        for edge in bg.get_edges_by_vertex(BGVertex(v)):
            # print(type(edge.vertex1), type(v))
            assert str(edge.vertex1) == v
            if get_colors_by_edge(edge)[BGGenome(genome)] == 1:
                return edge.vertex2

    def get_genome_color_by_edge(genome):
        if cnt[BGGenome(genome)] == 1:
            return 0
        else:
            v1, v2 = edge.vertex1.name, edge.vertex2.name
            if v1 > v2: v1, v2 = v2, v1

            # print(v1, v2)
            v1_neighbour = get_neighbour_with_genome(v1, genome)
            v2_neighbour = get_neighbour_with_genome(v2, genome)
            # print(v1_neighbour, v2_neighbour)
            if bg.get_edge_by_two_vertices(v1_neighbour, v2_neighbour):
                # print("hey!")
                pair = (v1_neighbour, v2_neighbour)
                if pair not in possible_edges:
                    possible_edges.append(pair)
                return 2 + possible_edges.index(pair)
            else:
                return 1

    cnt = get_colors_by_edge(edge)
    possible_edges = []
    return {genome: get_genome_color_by_edge(genome) for genome in genomes}, possible_edges


if __name__ == "__main__":
    labels_dict = make_labels_dict(csv_file)
    tree_holder = TreeHolder(tree_file, scale=scale, labels_dict=labels_dict,
                             colors=(
                                 'White', 'Gainsboro', 'LightGreen', 'LightBlue', 'NavajoWhite', 'LightPink',
                                 'LightCoral', 'Purple', 'Navy', 'Olive', 'Teal', 'SaddleBrown', 'SeaGreen', 'DarkCyan',
                                 'DarkOliveGreen', 'DarkSeaGreen'))

    genomes = get_genomes_contain_blocks(grimm_file)[0]
    print_species_stats(genomes, tree_holder.get_all_leafs())

    bg = GRIMMReader.get_breakpoint_graph(open(grimm_file))
    print('<bg parsed>')

    trees_folder = output_folder + 'trees/'
    os.makedirs(trees_folder, exist_ok=True)

    bg_folder = output_folder + 'bg_components/'
    os.makedirs(bg_folder, exist_ok=True)

    ans = [['vertex1', 'vertex2', 'parallel_rear_score', 'parallel_rear_unique_innovation_count',
            'parallel_rear_all_innovations_count', 'parallel_breakpoint_score', 'parallel_breakpoint_count']]

    # consistency_checker = TreeConsistencyChecker(tree_file)
    for i, component_bg in enumerate(bg.connected_components_subgraphs()):
        g = component_bg.bg
        nodes_len = len(list(component_bg.nodes()))
        if nodes_len == 2: continue

        print(f'component_{i}_size_{len(component_bg.bg)}')

        draw_bp_with_weights(component_bg, bg_folder + f'component_{i}_size_{len(component_bg.bg)}.pdf')
        # save_pygraphviz(component_bg, component_folder + 'graph_pygraphviz.pdf')

        for i_edge, edge in enumerate(component_bg.edges()):
            v1, v2 = edge.vertex1.name, edge.vertex2.name
            if v1 > v2: v1, v2 = v2, v1

            genome_colors, neighbour_edges = get_genome_colors_by_edge(component_bg, edge, genomes)

            # shigella differs?
            if white_proportion(genome_colors.values()) < 0.5: continue
            print(i_edge, 'of', len(list(component_bg.edges())), ':', v1, v2)

            # consistency
            tree_holder.count_innovations_fitch(genome_colors)

            score_rear, count_rear, count_all_rear = tree_holder.count_parallel_rearrangements(skip_grey=True)
            score_break, count_break = tree_holder.count_parallel_breakpoints()

            ans.append([v1, v2, score_rear, count_rear, count_all_rear, score_break, count_break])

            labels = ['edge exists', 'parallel edge doesn\'t exist'] + [f'parallel edge {v1}-{v2}' for (v1, v2) in
                                                                        neighbour_edges]

            tree_holder.draw(trees_folder + f'edge_{v1}_{v2}.pdf', legend_labels=labels, show_branch_support=True)

    with open(csv_out_file, 'w') as f:
        wtr = csv.writer(f)
        wtr.writerows(ans)