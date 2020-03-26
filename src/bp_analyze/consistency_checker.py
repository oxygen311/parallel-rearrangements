from Bio import Phylo
from collections import defaultdict, Counter

import networkx as nx

def find_root(G):
    for v in G.nodes():
        if G.in_degree(v) == 0:
            return v

def get_colorset(G,v, colorset):
    if v in colorset:
        return colorset[v]
    else:
        v1, v2 = G.neighbors(v)
        s1 = get_colorset(G,v1, colorset)
        s2 = get_colorset(G,v2, colorset)
        if len(s1 & s2) > 0:
            colorset[v] = s1 & s2
        else:
            colorset[v] = s1 | s2
        return colorset[v]


def down_to_leaves(G, v,  color, colorset, colors):
    if len(list(G.neighbors(v))) == 0:
        return 0
    for u in G.neighbors(v):
        if color in colorset[u]:
            colors[u] = color
        else:
            colors[u] = colorset[u].pop()
            colorset[u].add(colors[u])
        down_to_leaves(G,u,colors[u],colorset, colors)

def coloring(G,colors):
    #here only leaves are colored, colors is a dictionary of leaves colors
    colorset = defaultdict(set)
    for v in colors:
        colorset[v].add(colors[v])
    r = find_root(G)
    get_colorset(G,r,colorset)
    colors[r] = colorset[r].pop()
    down_to_leaves(G,r, colors[r] ,colorset,colors)

def is_convex(G, colors):
    #here only all the nodes are colored, colors is a dictionary of all the colors
    colored_nodes = defaultdict(list)
    for v in G.nodes():
        colored_nodes[colors[v]].append(v)

    res = True

    for color in colored_nodes:
        G_induced = G.subgraph(colored_nodes[color])
        if not(nx.is_weakly_connected(G_induced)):
            res = False
    return res

def recursive_convert(v, label, old_graph, new_graph):
    if len(list(old_graph.neighbors(v))) == 0: return

    v1, v2 = old_graph.neighbors(v)

    left_label = label + '0'  if v1.name is None else v1.name
    right_label = label + '1' if v2.name is None else v2.name

    new_graph.add_edge(label, left_label, length=v1.branch_length)
    new_graph.add_edge(label, right_label, length=v2.branch_length)

    recursive_convert(v1, left_label, old_graph, new_graph)
    recursive_convert(v2, right_label, old_graph, new_graph)

class TreeConsistencyChecker:
    def __init__(self, tree_file):
        phylo_tree = next(Phylo.parse(tree_file, 'newick', rooted=True))
        phylo_g = Phylo.to_networkx(phylo_tree)

        self.g = nx.DiGraph()
        root = find_root(phylo_g)

        recursive_convert(root, 's', phylo_g, self.g)

    def check_consistency(self, colors):
        # colors = colors_.copy()
        coloring(self.g, colors)
        return is_convex(self.g, colors)