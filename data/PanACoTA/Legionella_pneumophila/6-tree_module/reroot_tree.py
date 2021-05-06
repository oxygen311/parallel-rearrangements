from ete3 import Tree, TreeStyle

tree_file = 'LEPN_22.nucl.grp.aln.iqtree_tree.treefile'
t = Tree(tree_file)


circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
# circular_style.scale = 20

t.set_outgroup(t.get_common_ancestor(t & 'LEPN.0521.00010', t & 'LEPN.0521.00006'))
# t.set_outgroup(t.get_common_ancestor(t & 'STPY.0321.00009', t & 'STPY.0321.00002'))

t.render("test.pdf", tree_style=circular_style)
t.write(outfile=tree_file.replace('_tree', '_tree_reroot'))
