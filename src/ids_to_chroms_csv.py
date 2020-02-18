from glob import glob
from Bio import SeqIO

import os.path
import csv


folder = 'data/Bacillus_anthracis'
out_file = 'assemblies_chrs.csv'

ans = [['assembly', 'chromosome']]
for file in glob(f'{folder}/fna/*.fna'):
    if file.endswith('merged.fna'): continue

    contigs = [contig for contig in SeqIO.parse(open(file), 'fasta') if not 'plasmid' in contig.description]
    assert len(contigs) == 1

    ans.append([os.path.splitext(os.path.basename(file))[0], contigs[0].id])

with open(folder + '/' + out_file, 'w') as f:
    wr = csv.writer(f)
    wr.writerows(ans)