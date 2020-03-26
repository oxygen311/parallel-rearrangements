from glob import glob
from Bio import SeqIO

from collections import Counter

import os.path
import csv


folder = 'data/Staphylococcus_aureus'
out_file = 'assemblies_chrs.csv'

ans = [['assembly', 'chromosome', 'description']]
contigs_ids = []
unique_assembly_count = 0

for file in glob(f'{folder}/fna/*.fna'):
    if file.endswith('merged.fna'): continue
    unique_assembly_count += 1

    contigs = [contig for contig in SeqIO.parse(open(file), 'fasta') if not 'plasmid' in contig.description]
    assert len(contigs) == 1

    contigs_ids.append(contigs[0].id)

    ans.append([os.path.splitext(os.path.basename(file))[0],
                contigs[0].id,
                contigs[0].description.split(' ', 1)[1].replace(' chromosome, complete genome', '')
                                                       .replace(', complete genome', '')
                                                       .replace(', complete sequence', '')
                                                       .replace(' genome assembly, chromosome: 1', '')
                                                       .replace(' complete genome', '')])

with open(folder + '/' + out_file, 'w') as f:
    wr = csv.writer(f)
    wr.writerows(ans)

print("Unique assemblies:", unique_assembly_count)

cnt_contigs = Counter(contigs_ids)

for contig, count in cnt_contigs.items():
    if count > 1:
        print(f'WARNING!!! contig {contig} appears {count} times')