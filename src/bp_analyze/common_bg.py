from bg.breakpoint_graph import BreakpointGraph

import matplotlib.pyplot as plt
import networkx as nx

get_colors_by_edge = lambda e: e.multicolor.multicolors
get_weight_by_edge = lambda e: sum(get_colors_by_edge(e).values()) / 15

def get_neighbors_subgraph(bg, edge):
    g = BreakpointGraph()
    [g.add_bgedge(n_edge) for n_edge in bg.get_edges_by_vertex(edge.vertex1) if n_edge != edge]
    [g.add_bgedge(n_edge) for n_edge in bg.get_edges_by_vertex(edge.vertex2) if n_edge != edge]
    g.add_bgedge(edge)
    return g

def draw_bp_with_weights(g, tree, file=None, color='red', do_show=False):
    def get_color(edge):
        if not tree.bgedge_is_tree_consistent(edge): return color
        return 'black'

    g_nx = nx.Graph()
    [g_nx.add_edge(edge.vertex1.name, edge.vertex2.name, weight=get_weight_by_edge(edge), color=get_color(edge)) for
     edge in g.edges()]

    pos = nx.circular_layout(g_nx)

    edges = g_nx.edges()
    colors = [g_nx[u][v]['color'] for u, v in edges]
    weights = [g_nx[u][v]['weight'] for u, v in edges]

    plt.figure()
    nx.draw(g_nx, pos, edges=edges, edge_color=colors, width=weights, with_labels=True, node_size=500,
            node_color='#B6D3E6')
    if not file is None: plt.savefig(file)
    if do_show: plt.show()