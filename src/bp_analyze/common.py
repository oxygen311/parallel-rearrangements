from bg.grimm import GRIMMReader
from collections import defaultdict, Counter

from operator import itemgetter

import pandas as pd

def make_labels_dict(file, row_from='chromosome', row_to='description'):
    try:
        df = pd.read_csv(file)
        return {row[row_from]: row[row_to] for i, row in df.iterrows()}
    except FileNotFoundError:
        return None

def get_genomes_contain_blocks(grimm_file, skip_unique=False, split_by=None):
    genomes, blocks = set(), set()

    with open(grimm_file) as f: ls = f.readlines()
    block_genome_count = defaultdict(Counter)

    for i in range(0, len(ls), 2):
        name = GRIMMReader.parse_genome_declaration_string(ls[i]).name
        if split_by is not None:
            name = name.split(split_by)[0]
        data = GRIMMReader.parse_data_string(ls[i + 1])[1]
        genomes.add(name)
        for _, block in data:
            blocks.add(int(block))
            block_genome_count[int(block)][name] += 1

    print('Got blocks:', len(blocks))
    if skip_unique:

        blocks = [block for block in blocks
                  if not all(block_genome_count[block][genome] == 1 for genome in genomes)]
        print('After skipping unique:', len(blocks))

    return genomes, list(sorted(blocks)), block_genome_count


def get_block_seqs(grimm_file):
    with open(grimm_file) as f: ls = f.readlines()

    seqs = []
    for i in range(0, len(ls), 2):
        data = list(map(itemgetter(1), GRIMMReader.parse_data_string(ls[i + 1])[1]))
        seqs.append(data)

    return seqs

def print_species_stats(block_genomes, tree_genomes):
    print('Block species:', len(block_genomes))
    print('Tree species:', len(tree_genomes))
    print('Intersected species:', len(block_genomes & tree_genomes))
    print()

    print(block_genomes - (block_genomes & tree_genomes))
    print(tree_genomes - (block_genomes & tree_genomes))
    print()

# def group_by_count(genomes_count):
#     ans = defaultdict(set)
#     for genome, count in genomes_count.items():
#         ans[count].add(genome)
#
#     return ans.values()