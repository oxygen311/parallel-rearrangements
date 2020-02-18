from bg.grimm import GRIMMReader
from bg.tree import BGTree
from Bio import Phylo

import os
import pylab

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

csv_data_file = 'data/12862_2019_1403_MOESM1_ESM.csv'

color = 'purple'
sp_folder = 'sibeliaz/Streptococcus_pyogenes/'
blocks_folder = sp_folder + 'fine/5000/'
grimm_file = blocks_folder + 'genomes_permutations_unique.txt'
tree_file = sp_folder + 'RAxML_bestTree.concat_alignment.newick'
tree_support_file = sp_folder + 'RAxML_support.newick'
output_folder = blocks_folder + 'bp_graph_output'

get_colors_by_edge = lambda e: e.multicolor.multicolors
get_weight_by_edge = lambda e: sum(get_colors_by_edge(e).values()) / 7


def construct_color_func(cs):
    def func(c):
        if c in cs: return color
        return 'black'

    return func


def construct_label_func():
    df = pd.read_csv('data/12862_2019_1403_MOESM1_ESM.csv')
    strains = {row['ACCESSION']: row['Strain'] for index, row in df.iterrows()}
    return lambda id: strains.get(str(id), '')


label_func = construct_label_func()


def draw_bp_with_weights(g, tree, file=None):
    def get_color(edge):
        if not tree.bgedge_is_tree_consistent(edge): return color
        # if not tree.bgedge_is_vtree_consistent(edge): return 'blue'
        return 'black'

    g_nx = nx.Graph()
    [g_nx.add_edge(edge.vertex1.name, edge.vertex2.name, weight=get_weight_by_edge(edge), color=get_color(edge)) for
     edge in g.edges()]

    pos = nx.circular_layout(g_nx)

    edges = g_nx.edges()
    colors = [g_nx[u][v]['color'] for u, v in edges]
    weights = [g_nx[u][v]['weight'] for u, v in edges]

    nx.draw(g_nx, pos, edges=edges, edge_color=colors, width=weights, with_labels=True, node_size=500,
            node_color='#B6D3E6')
    if not file is None: plt.savefig(file)
    plt.show()


bg = GRIMMReader.get_breakpoint_graph(open(grimm_file))
tree_string = open(tree_file).read()

phylo_tree = next(Phylo.parse(tree_file, 'newick'))
phylo_support_tree = next(Phylo.parse(tree_support_file, 'newick'))

tree = BGTree(tree_string)
for i, component_bg in enumerate(bg.connected_components_subgraphs()):
    g = component_bg.bg
    nodes_len = len(list(component_bg.nodes()))
    if nodes_len == 2: continue

    # draw graph
    component_folder = output_folder + f'component_{i}/'
    if not os.path.exists(component_folder): os.makedirs(component_folder)
    draw_bp_with_weights(component_bg, tree, component_folder + 'graph.pdf')

    # draw not tree consistent edges
    for edge in component_bg.edges():
        if not tree.bgedge_is_tree_consistent(edge):
            plt.rcParams.update({'font.size': 6, 'lines.linewidth': 0.4})
            species = [label_func(color.name) for color in get_colors_by_edge(edge).keys()]
            Phylo.draw(phylo_support_tree, label_colors=construct_color_func(species), do_show=False,
                       label_func=label_func)
            # plt.rcParams.update({'font.size': 14, 'lines.linewidth': 1})
            # species = [color.name for color in get_colors_by_edge(edge).keys()]
            # Phylo.draw(phylo_tree, label_colors=construct_color_func(species), do_show=False)
            plt.axis('off')

            v1, v2 = edge.vertex1.name, edge.vertex2.name
            if v1 > v2: v1, v2 = v2, v1
            plt.savefig(component_folder + f'edge_{v1}_{v2}.pdf', bbox_inches='tight')
            print("Saved to ", component_folder + f'edge_{v1}_{v2}.pdf')
            plt.show()
