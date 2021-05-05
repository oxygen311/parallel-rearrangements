from ete3 import Tree, TreeStyle

t = Tree('asm_tree.tree')


circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
# circular_style.scale = 20

# t.set_outgroup(t.get_common_ancestor(t & 'GCA_002854065.1', t & 'GCA_000007445.1'))
t.set_outgroup(t.get_common_ancestor(t & 'GCA_002310615.1', t & 'GCA_000007445.1'))

t.write(outfile='asm_tree_reroot.tree')
# t.render("test.pdf", tree_style=circular_style)
