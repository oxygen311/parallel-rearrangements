import os

import pandas as pd

from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks, get_block_seqs
from src.bp_analyze.tree_holder import TreeHolder

block_size = 5000
folder = f'data/EPEC/blocks/{block_size}/'
grimm_file = 'genomes_permutations.txt'
tree_file = 'data/EPEC/RAxML_bestTree.whole_tree'

csv_file = 'data/EPEC/aEPECmeta.csv'

row_of_interest = 'inflammatory_potential_in_vitro'
# row_of_interest = 'disease_healthy'

gain_csv = f'information_gain_{row_of_interest}.csv'

best_count = 42
output_folder = folder + f'best_{best_count}_trees_{row_of_interest}/'

# type_colors = {
#     'healthy': 'Black',
#     'disease': 'Crimson',
#     'unknown': 'Gray'
# }

type_colors = {
    'nan': 'Black',
    'pro_inflammatory': 'Crimson',
    'anti_inflammatory': 'Purple',
}

def make_node_colors():
    df_types = pd.read_csv(csv_file)
    return {'EPEC' + str(row.EPEC_number): type_colors[str(row[row_of_interest])] for _, row in df_types.iterrows()}


if __name__ == "__main__":
    genomes, blocks, block_genome_count = get_genomes_contain_blocks(folder + grimm_file, skip_unique=True,
                                                                     split_by='.')

    node_colors = make_node_colors()
    tree_holder = TreeHolder(tree_file, scale=4200, node_colors=node_colors)

    df_gain = pd.read_csv(folder + gain_csv)
    df_gain.sort_values('gain', ascending=False, inplace=True)

    os.makedirs(output_folder, exist_ok=True)

    for i, row in zip(range(best_count), df_gain.itertuples()):
        count = block_genome_count[row.block]
        tree_holder.count_innovations_fitch({genome: count[genome] for genome in genomes})

        labels = [f'{i}{"+" if i == tree_holder.max_color - 1 else ""} copies'
                  for i in range(max(count.values()) + 1)]

        tree_holder.draw(output_folder + f'{i}_block_{row.block}_gain_{"{:.4f}".format(row.gain)}.pdf',
                         legend_labels=labels, show_branch_support=False, mode='r', legend_scale=0.3,
                         node_colors=type_colors)
