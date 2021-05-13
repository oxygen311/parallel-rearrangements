from ete3 import Tree, TreeStyle

from glob import glob

tree_file = glob('*.treefile')[0]
t = Tree(tree_file)


circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
# circular_style.scale = 20

# t.set_outgroup(t.get_common_ancestor(t & 'LEPN.0521.00010', t & 'LEPN.0521.00006'))
t.set_outgroup(t.get_common_ancestor(t & 'LIMO.0521.00128', t & 'LIMO.0521.00247'))

t.prune(set(t.get_leaves()) - set([t & 'LIMO.0521.00221']))

t.render("test2.pdf", tree_style=circular_style)
t.write(outfile=tree_file.replace('_tree', '_tree_reroot'))
