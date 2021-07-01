from bg.grimm import GRIMMReader
from bg.vertices import BGVertex
from bg.tree import BGTree

from src.bp_analyze.tree_holder import TreeHolder
from src.bp_analyze.common import get_genomes_contain_blocks, print_species_stats, make_labels_dict
from src.bp_analyze.call_unique_gene import get_genome_colors_by_edge

from pprint import pprint

sp_folder = 'data/E_coli/'
testing = False
csv_file = sp_folder + 'total_stats.csv'
scale = 15000
blocks_folder = sp_folder + 'sibeliaz_out/fine/1000/'
grimm_file = blocks_folder + 'genomes_permutations_unique.txt'
tree_file = sp_folder + 'tree/tree_midpoint_converted.nwk'
output_folder = blocks_folder + 'bp_unique_block_output/'


if __name__ == "__main__":
    labels_dict = make_labels_dict(csv_file)
    tree_drawer = TreeHolder(tree_file, scale=scale, labels_dict=labels_dict,
                             colors=(
                                 'White', 'Gainsboro', 'LightGreen', 'LightBlue', 'NavajoWhite', 'LightPink',
                                 'LightCoral', 'Purple', 'Navy', 'Olive', 'Teal', 'SaddleBrown', 'SeaGreen', 'DarkCyan',
                                 'DarkOliveGreen', 'DarkSeaGreen'))

    genomes = get_genomes_contain_blocks(grimm_file)[0]

    bg = GRIMMReader.get_breakpoint_graph(open(grimm_file))
    print('<bg parsed>')

    v1, v2 = '475h', '477t'

    edge = bg.get_edge_by_two_vertices(BGVertex(v1), BGVertex(v2))
    print(len(list(bg.edges())))
    # print(bg, type(bg))
    # print(edge, type(edge))
    # print(genomes, type(genomes))
    genome_colors, neighbour_edges = get_genome_colors_by_edge(bg, edge, genomes)

    pprint(genome_colors)
    labels = ['edge exists', 'parallel edge doesn\'t exist'] + [f'Inversion {v1}-{v2}-{v3}-{v4}' for (v3, v4) in neighbour_edges]
    #
    tree_drawer.draw(f'edge_{v1}_{v2}.pdf', colors_dict=genome_colors, show_branch_support=False,
                     legend_labels=labels, show_scale=False, legend_scale=1, draw_pres=True)
    # tree_drawer.draw(f'tree.pdf', colors_dict=defaultdict(int))