from bg.grimm import GRIMMReader
from bg.tree import BGTree
from bg.genome import BGGenome
from bg.vertices import BGVertex

get_colors_by_edge = lambda e: e.multicolor.multicolors

bg = GRIMMReader.get_breakpoint_graph(open('genomes_permutations_unique.txt'))
for v in bg.bg:
    edges_counts = [len(get_colors_by_edge(edge)) for edge in bg.get_edges_by_vertex(v)]

    print(v, edges_counts)
    # print(list(map(lambda x: x > sum(edges_counts) / 2, edges_counts)))

    assert any(map(lambda x: x > sum(edges_counts) / 2, edges_counts))
