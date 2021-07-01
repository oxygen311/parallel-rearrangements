from src.bp_analyze.tree_holder import TreeHolder

from ete3 import Tree

tree = '(((Strain 1,Strain 2),Strain 3:2),(Strain 4,Strain 5));'
# tree = '((Strain 1:0.1,Strain 2:0.1):0.2,((Strain 3:0.1, Strain 4:0.1):0.1,Strain 5:0.2):0.1);'

colors = {
    'Strain 1': 0,
    'Strain 2': 1,
    'Strain 3': 1,
    'Strain 4': 2,
    'Strain 5': 3,
}

node_colors = {
    'Strain 1': '#0070ba',
    'Strain 2': '#fe2e59',
    'Strain 3': '#813891',
    'Strain 4': '#239e5c',
    'Strain 5': '#777777'
}

out_filename = 'tree_example_3.pdf'

td = TreeHolder(tree, scale=42, node_colors=node_colors, colors=(
                                 'White', '#ffdddd', '#cbf7cb', '#eeeeee', 'NavajoWhite', 'LightPink'))

labels = ['Edge exists', 'Inversion 1h-2t-4h-5t', 'Inversion 1h-2t-5h-0t', 'Complex break']

td.count_innovations_fitch(colors)
td.draw(out_filename, show_branch_support=False, color_internal_nodes=False, show_scale=False, legend_labels=labels, legend_scale=0.25, draw_pres=True)