from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace

class TreeDrawer():

    def __init__(self, tree, scale, labels_dict=None, colors=('Gainsboro', 'White', 'LightGreen', 'LightBlue', 'LightPink', 'LightCoral')):
        self.tree = Tree(tree)
        self.scale = scale
        self.colors = colors
        self.max_color = len(self.colors)

        for node in self.tree.traverse():
            # Hide node circles
            node.img_style['size'] = 0

            if node.is_leaf():
                name_face = TextFace(labels_dict[node.name] if labels_dict else node.name, fsize=10)
                node.add_face(name_face, column=0)

    def draw(self, file, colors_dict, color_internal_nodes=True, legend_labels=[]):
        # tree = self.tree.copy()
        for node in self.tree.traverse():
            if node.is_leaf():
                node.img_style['bgcolor'] = self.colors[min(colors_dict[node.name], self.max_color - 1)]
            elif color_internal_nodes:
                leaf_colors = set(colors_dict[leaf] for leaf in node.iter_leaf_names())
                if len(leaf_colors) == 1:
                    node.img_style['bgcolor'] = self.colors[min(leaf_colors.pop(), self.max_color - 1)]
                else:
                    node.img_style['bgcolor'] = 'White'

        # Set up style for circular tree
        ts = TreeStyle()
        ts.mode = "c"
        ts.scale = self.scale
        # Disable the default tip names config
        ts.show_leaf_name = False
        ts.show_branch_support = True

        current_colors = self.colors[0:max(colors_dict.values()) + 1]
        # if annotate == 'dups':
        #     legend_labels = [f'{i}{"+" if i == len(colors) - 1 else ""} copies' for i in range(len(current_colors))]
        # elif annotate == 'unique':
        #     legend_labels =
        # else:
        #     legend_labels = []

        for i, (label, color) in enumerate(zip(legend_labels, current_colors)):
            ts.legend.add_face(CircleFace(18, color), column=0)
            ts.legend.add_face(CircleFace(10, 'White'), column=1)
            ts.legend.add_face(TextFace(label, fsize=40), column=2)
            ts.legend.add_face(CircleFace(10, 'White'), column=3)

        # self.tree.render("ete_tree.pdf", dpi=300, tree_style=ts)
        self.tree.render(file, w=1000, tree_style=ts)

    def get_all_leafs(self):
        return {node.name for node in self.tree.traverse() if node.is_leaf()}