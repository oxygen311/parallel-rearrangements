from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace

from collections import defaultdict


class TreeDrawer():

    def __init__(self, tree, scale, labels_dict=None, colors=(
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

    def draw(self, file, colors_dict, color_internal_nodes=True, legend_labels=(), show_branch_support=True,
             show_scale=True, legend_scale=1, draw_pres=False):
        # tree = self.tree.copy()
        for node in self.tree.traverse():
            if node.is_leaf():
                color = self.colors[min(colors_dict[node.name], self.max_color - 1)]
                if color != 'White' or not draw_pres:
                    node.img_style['bgcolor'] = color

            elif color_internal_nodes:
                leaf_colors = set(colors_dict[leaf] for leaf in node.iter_leaf_names())
                if len(leaf_colors) == 1:
                    color = self.colors[min(leaf_colors.pop(), self.max_color - 1)]
                    if color != 'White' or not draw_pres:
                        node.img_style['bgcolor'] = color

        # Set up style for circular tree
        ts = TreeStyle()
        ts.mode = "c"
        ts.scale = self.scale
        # Disable the default tip names config
        ts.show_leaf_name = False
        # ts.force_topology = True
        ts.show_branch_support = show_branch_support

        # ts.branch_vertical_margin = 20
        ts.show_scale = show_scale
        current_colors = self.colors[0:max(colors_dict.values()) + 1]

        for i, (label, color_) in enumerate(zip(legend_labels, current_colors)):
            ts.legend.add_face(CircleFace(24 * legend_scale, color_), column=0)
            ts.legend.add_face(CircleFace(13 * legend_scale, 'White'), column=1)
            ts.legend.add_face(TextFace(label, fsize=53 * legend_scale), column=2)
            ts.legend.add_face(CircleFace(13 * legend_scale, 'White'), column=3)

        # self.tree.render("ete_tree.pdf", dpi=300, tree_style=ts)
        self.tree.render(file, w=1000, tree_style=ts)

    def get_all_leafs(self):
        return {node.name for node in self.tree.traverse() if node.is_leaf()}