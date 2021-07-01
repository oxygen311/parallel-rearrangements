from bg.grimm import GRIMMReader
from bg.tree import BGTree
from bg.genome import BGGenome
from bg.vertices import BGVertex

from src.bp_analyze.tree_holder import TreeHolder
from src.bp_analyze.common import get_genomes_contain_blocks, print_species_stats, make_labels_dict
from src.bp_analyze.common_bg import draw_bp_with_weights, get_colors_by_edge, save_pygraphviz
from src.bp_analyze.consistency_checker import TreeConsistencyChecker
from src.bp_analyze.shigella_common import get_class_colors, white_proportion, count_shigella_differs_from_white

from pprint import pprint

import os

sp_folder = 'data/E_coli/'
testing = False
csv_file = sp_folder + 'total_stats.csv'
scale = 15000
blocks_folder = sp_folder + 'sibeliaz_out/fine/1000/'
grimm_file = blocks_folder + 'genomes_permutations_unique.txt'
tree_file = sp_folder + 'tree/tree_midpoint_converted.nwk'
output_folder = blocks_folder + 'bp_unique_block_output/'


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
    tree_drawer = TreeHolder(tree_file, scale=scale, labels_dict=labels_dict,
                             colors=(
                                 'White', 'Gainsboro', 'LightGreen', 'LightBlue', 'NavajoWhite', 'LightPink',
                                 'LightCoral', 'Purple', 'Navy', 'Olive', 'Teal', 'SaddleBrown', 'SeaGreen', 'DarkCyan',
                                 'DarkOliveGreen', 'DarkSeaGreen'))

    genomes = get_genomes_contain_blocks(grimm_file)[0]
    print_species_stats(genomes, tree_drawer.get_all_leafs())

    bg = GRIMMReader.get_breakpoint_graph(open(grimm_file))
    bg_tree = BGTree(tree_file)
    print('<bg parsed>')

    consistency_checker = TreeConsistencyChecker(tree_file)
    for i, component_bg in enumerate(bg.connected_components_subgraphs()):
        g = component_bg.bg
        nodes_len = len(list(component_bg.nodes()))
        if nodes_len == 2: continue

        # draw graph
        component_folder = output_folder + f'component_{i}_size_{len(component_bg.bg)}/'
        print(f'component_{i}_size_{len(component_bg.bg)}')
        os.makedirs(component_folder, exist_ok=True)

        draw_bp_with_weights(component_bg, bg_tree, component_folder + 'graph_networkx.pdf')
        # save_pygraphviz(component_bg, component_folder + 'graph_pygraphviz.pdf')

        for i_edge, edge in enumerate(component_bg.edges()):
            v1, v2 = edge.vertex1.name, edge.vertex2.name
            if v1 > v2: v1, v2 = v2, v1

            print(i_edge, 'of', len(list(component_bg.edges())), ':', v1, v2)
            genome_colors, neighbour_edges = get_genome_colors_by_edge(component_bg, edge, genomes)

            # shigella differs?
            class_colors = get_class_colors(labels_dict, genome_colors)
            if white_proportion(class_colors['Escherichia coli']) < 0.5: continue
            shigella_differs = count_shigella_differs_from_white(class_colors)

            # consistency
            type = 'consistent'
            inconsistent_colors = consistency_checker.check_consistency(genome_colors.copy())
            if len(inconsistent_colors) == 0:
                if len(consistency_checker.check_consistency({k: (0 if v == 0 else 1) for k, v in genome_colors.items()})) > 0:
                    type = 'parallel_breakpoint'
            else:
                type = 'parallel_' + \
                ('breakpoint' if len(inconsistent_colors) == 1 and inconsistent_colors.pop() == 1 else 'rearrangement')

            curr_output_folder = component_folder + f'{type}__{shigella_differs}_shigells_differ/'
            os.makedirs(curr_output_folder, exist_ok=True)

            labels = ['edge exists', 'parallel edge doesn\'t exist'] + [f'parallel edge {v1}-{v2}' for (v1, v2) in
                                                                        neighbour_edges]

            tree_drawer.draw(curr_output_folder + f'edge_{v1}_{v2}.pdf', colors_dict=genome_colors,
                             legend_labels=labels)