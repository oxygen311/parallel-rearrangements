from bg.grimm import GRIMMReader
from collections import defaultdict, Counter


import pandas as pd

def make_labels_dict(file):
    try:
        df = pd.read_csv(file)
        return {row['AssemblyID']: row['Organism'] for i, row in df.iterrows()}
    except FileNotFoundError:
        return None

def get_genomes_contain_blocks(grimm_file):
    genomes, blocks = set(), set()

    with open(grimm_file) as f: ls = f.readlines()
    block_genome_count = defaultdict(Counter)

    for i in range(0, len(ls), 2):
        name = GRIMMReader.parse_genome_declaration_string(ls[i]).name
        data = GRIMMReader.parse_data_string(ls[i + 1])[1]
        genomes.add(name)
        for _, block in data:
            blocks.add(block)
            block_genome_count[block][name] += 1

    return genomes, list(blocks), block_genome_count

def print_species_stats(block_genomes, tree_genomes):
    print('Block species:', len(block_genomes))
    print('Tree species:', len(tree_genomes))
    print('Intersected species:', len(block_genomes & tree_genomes))
    print()

# def group_by_count(genomes_count):
#     ans = defaultdict(set)
#     for genome, count in genomes_count.items():
#         ans[count].add(genome)
#
#     return ans.values()