from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace, RectFace

from collections import defaultdict, Counter
from itertools import combinations

from src.bp_analyze.common import get_genomes_contain_blocks, print_species_stats, make_labels_dict


class TreeHolder:
    def __init__(self, tree, scale=None, labels_dict=None, colors=(
            'Gainsboro', 'White', 'LightGreen', 'LightBlue', 'NavajoWhite', 'LightPink',
            'LightCoral', 'Purple', 'Navy', 'Olive', 'Teal', 'SaddleBrown', 'SeaGreen', 'DarkCyan',
            'DarkOliveGreen', 'DarkSeaGreen'), node_colors=defaultdict(lambda: 'black')):
        self.tree = Tree(tree)
        self.scale = scale
        self.colors = colors
        self.max_color = len(self.colors)

        for node in self.tree.traverse():
            # Hide node circles
            node.img_style['size'] = 0

            if node.is_leaf():
                name_face = TextFace(labels_dict[node.name] if labels_dict else node.name,
                                     fgcolor=node_colors[node.name])
                node.add_face(name_face, column=0)

    def draw(self, file, color_internal_nodes=True, legend_labels=(), show_branch_support=True,
             show_scale=True, legend_scale=1, draw_pres=False, mode="c", node_colors={}):
        # tree = self.tree.copy()
        for node in self.tree.traverse():
            if not (color_internal_nodes or node.is_leaf()): continue
            color = self.colors[min(node.color, self.max_color - 1)]
            if not (color == 'White' and draw_pres):
                node.img_style['bgcolor'] = color

        ts = TreeStyle()
        ts.mode = mode
        ts.scale = self.scale
        # Disable the default tip names config
        ts.show_leaf_name = False
        ts.show_branch_support = show_branch_support

        # ts.branch_vertical_margin = 20
        ts.show_scale = show_scale
        cur_max_color = max(v.color for v in self.tree.traverse())
        current_colors = self.colors[0:cur_max_color + 1]

        for label, color_ in zip(legend_labels, current_colors):
            ts.legend.add_face(CircleFace(24 * legend_scale, color_), column=0)
            ts.legend.add_face(CircleFace(13 * legend_scale, 'White'), column=1)
            ts.legend.add_face(TextFace(label, fsize=53 * legend_scale), column=2)
            ts.legend.add_face(CircleFace(13 * legend_scale, 'White'), column=3)

        for node, color_ in node_colors.items():
            ts.legend.add_face(CircleFace(24 * legend_scale, 'White'), column=4)
            ts.legend.add_face(TextFace(node, fsize=53 * legend_scale, fgcolor=color_), column=5)

        # self.tree.render("ete_tree.pdf", dpi=300, tree_style=ts)
        self.tree.render(file, w=1000, tree_style=ts)

    def get_all_leafs(self):
        return {node.name for node in self.tree.traverse() if node.is_leaf()}

    def count_innovations_fitch(self, leaf_colors):
        def assign_colorset_feature(v):
            if v.is_leaf():
                v.add_features(colorset={leaf_colors[v.name]}, color=leaf_colors[v.name])
            else:
                try:
                    child1, child2 = v.children
                except ValueError:
                    print(v.children)
                    raise ValueError('Tree must me binary')
                cs1 = assign_colorset_feature(child1)
                cs2 = assign_colorset_feature(child2)
                v.add_features(colorset=(cs1 & cs2) if len(cs1 & cs2) > 0 else cs1 | cs2)

            return v.colorset

        def chose_color(colorset):
            return sorted(colorset, key=lambda c: color_counter[c], reverse=True)[0]

        def down_to_leaves(v, color):
            if v.is_leaf(): return
            v.add_features(color=color if color in v.colorset else chose_color(v.colorset))
            for child in v.children:
                down_to_leaves(child, v.color)

        def count_innovations(v, innovations):
            for child in v.children:
                if v.color != child.color:
                    innovations[child.color].append(child)
                count_innovations(child, innovations)

        color_counter = Counter(leaf_colors.values())

        # get colorsets for internal nodes
        root = self.tree.get_tree_root()
        assign_colorset_feature(root)

        # get color for internal nodes
        root_color = chose_color(root.colorset)
        down_to_leaves(root, root_color)

        # get inconsistent colors
        self.innovations = defaultdict(list)
        count_innovations(root, self.innovations)

    def count_parallel_rearrangements(self, skip_grey):
        score, count, count_all = 0, 0, 0
        for color, nodes in self.innovations.items():
            if len(nodes) <= 1 or (skip_grey and color == 1): continue
            count += 1
            count_all += len(nodes)
            for n1, n2 in combinations(nodes, 2):
                score += n1.get_distance(n2)
        return score, count, count_all

    def count_parallel_breakpoints(self):
        count = sum(map(len, self.innovations.values()))
        score = sum(
            n1.get_distance(n2) for n1, n2 in combinations((n for ns in self.innovations.values() for n in ns), 2))
        return score, count

    def draw_coloring(self, file):
        for node in self.tree.traverse():
            node.img_style['bgcolor'] = self.colors[node.color]
        ts = TreeStyle()
        ts.show_leaf_name = False
        self.tree.render(file, w=1000, tree_style=ts)


# csv_file = 'data/Streptococcus_pyogenes/assemblies_chrs.csv'
# labels_dict = make_labels_dict(csv_file)
#
# t = 'data/Staphylococcus_aureus/tree/RAxML_bipartitionsBranchLabels.attepmt_converted_renamed'
# # t = 'data/Streptococcus_pyogenes/tree/RAxML_bestTree.concat_alignment.tree'
# th = TreeHolder(t, colors=('#FAFAFA', 'LightGreen'), labels_dict=None)
#
# d = defaultdict(int)
# # d = defaultdict(lambda: 0, {'CP007537': 1, 'CP014027': 1, 'CP008776': 1, 'AP012491': 1, 'AM295007': 1, 'CP008695': 1,
# #                             'CP007562': 1, 'AE014074': 1, 'BA000034': 1, 'CP007561': 1, 'CP010449': 1, 'CP014139': 1,
# #                             'CP011535': 1})
# #
# th.count_innovations_fitch({genome: d[genome] for genome in th.get_all_leafs()})
#
# th.draw('test.pdf', show_branch_support=True, color_internal_nodes=True, draw_pres=False)

csv_file = 'data/Staphylococcus_aureus/assemblies_chrs.csv'
labels_dict = make_labels_dict(csv_file)

t = 'data/Staphylococcus_aureus/tree/RAxML_bipartitionsBranchLabels.attepmt_converted_renamed'
t_2 = 'data/Staphylococcus_aureus/tree/RAxML_bipartitionsBranchLabels.attepmt_converted_renamed_midpoint'
# t = 'data/Streptococcus_pyogenes/tree/RAxML_bestTree.concat_alignment.tree'
th = TreeHolder(t, colors=('White', 'LightGreen'), labels_dict=labels_dict, scale=72000)

d = defaultdict(int)

# node = th.tree.get_common_ancestor('NC_017338.1', 'NZ_CP038460.1')
# th.tree.set_outgroup(node)
# th.tree.write(outfile=t_2)

th.count_innovations_fitch({genome: d[genome] for genome in th.get_all_leafs()})
th.draw('test.pdf', show_branch_support=True, color_internal_nodes=True, draw_pres=False)