from bg.grimm import GRIMMReader
from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx

grimm_file = 'data/Legionella_pneumophila/sibeliaz/fine/5000/genomes_permutations.txt'
# grimm_file = 'fine_sibelia_output/Streptococcus_pneumoniae/genomes_permutations.txt'

lines = open(grimm_file).read().split('\n')

class Unique_Filter:
    def __init__(self):
        self.allowed_blocks = {}
        self.first_call = True

    def update_allowed_blocks(self, ps):
        vs = [p[1] for p in ps]
        if self.first_call == True:
            self.allowed_blocks = vs
            self.first_call = False
        counter = Counter(vs)
        # print({b for b in vs if counter[b] > 2})
        self.allowed_blocks = {b for b in self.allowed_blocks if counter[b] == 1}

    def filter_unique(self, ps):
        return [p for p in ps if p[1] in self.allowed_blocks]

# make unique blocks list
i = 0
flt = Unique_Filter()
while i < len(lines):
    line = lines[i]
    if GRIMMReader.is_genome_declaration_string(line):
        data_line = lines[i + 1]
        parsed = GRIMMReader.parse_data_string(data_line)[1]
        flt.update_allowed_blocks(parsed)
        i += 2
    else:
        i += 1

print('Number of unique blocks:', len(flt.allowed_blocks))
# write allowed blocks
i = 0
with open(grimm_file.replace('.txt', '_unique.txt'), 'w') as f:
    while i < len(lines):
        line = lines[i]
        if GRIMMReader.is_genome_declaration_string(line):
            data_line = lines[i + 1]

            parsed = GRIMMReader.parse_data_string(data_line)[1]
            parsed = flt.filter_unique(parsed)

            print(line.split('.')[0], file=f)
            print(' '.join(p[0] + p[1] for p in parsed), '@', file=f)
            i += 2
        else:
            i += 1
