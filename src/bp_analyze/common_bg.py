from bg.breakpoint_graph import BreakpointGraph

import matplotlib.pyplot as plt
import networkx as nx
import pygraphviz as pgv

import math


get_colors_by_edge = lambda e: e.multicolor.multicolors
get_weight_by_edge = lambda e: sum(get_colors_by_edge(e).values()) / 15

def get_neighbors_subgraph(bg, edge):
    g = BreakpointGraph()
    [g.add_bgedge(n_edge) for n_edge in bg.get_edges_by_vertex(edge.vertex1) if n_edge != edge]
    [g.add_bgedge(n_edge) for n_edge in bg.get_edges_by_vertex(edge.vertex2) if n_edge != edge]
    g.add_bgedge(edge)
    return g

def draw_bp_with_weights(g, file=None, color='red', do_show=False):
    def get_color(edge):
        # if not tree.bgedge_is_tree_consistent(edge): return color
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

def save_pygraphviz(g, filename, prog='neato', max_parallel_edges=14):
    ind_by_node = lambda n: int(n[:-1]) * 2 + (1 if n[-1] == 'h' else 0)

    def cacl_pos(node):
        radius = n ** (6 / 7) / 2
        ind = sorted_ns.index(node)
        return "%f,%f!" % (math.cos(math.pi / n * ind * 2) * radius,
                           math.sin(math.pi / n * ind * 2) * radius)

    a = pgv.AGraph(strict=False)
    n = len(list(g.nodes()))

    colors = ['aquamarine', 'bisque', 'blue', 'blueviolet', 'brown', 'cadetblue1', 'chartreuse', 'chocolate', 'coral',
              'darkgoldenrod1', 'darkkhaki', 'darkolivegreen', 'darkorchid', 'darkslategray', 'deeppink', 'firebrick1',
              'forestgreen', 'goldenrod']

    ns = list(map(lambda c: (int(c.name[:-1]), 'a' if c.name[-1] == 't' else 'b'), g.nodes()))
    sorted_ns = [str(n) + ('t' if o == 'a' else 'h') for n, o in sorted(ns)]

    for node in g.nodes():
        # print(node)
        a.add_node(str(node), shape="circle", pos=cacl_pos(str(node)))

    for edge in g.edges():
        for sp, _ in zip(get_colors_by_edge(edge).keys(), range(max_parallel_edges)):
            a.add_edge(edge.vertex1.name, edge.vertex2.name, color=colors[hash(sp) % len(colors)])

    a.draw(filename, prog=prog)