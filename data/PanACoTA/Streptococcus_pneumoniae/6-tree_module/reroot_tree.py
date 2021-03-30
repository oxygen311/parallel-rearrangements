from ete3 import Tree, TreeStyle


t = Tree('STPN_97.nucl.grp.aln.iqtree_tree.treefile')


circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
# circular_style.scale = 20

# t.set_outgroup(t.get_common_ancestor(t & 'GCA_002854065.1', t & 'GCA_000007445.1'))
t.set_outgroup(t.get_common_ancestor(t & 'STPN.0321.00097', t & 'STPN.0321.00002'))

t.render("test.pdf", tree_style=circular_style)
t.write(outfile='STPN_97.nucl.grp.aln.iqtree_tree_reroot.treefile')
