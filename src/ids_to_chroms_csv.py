from glob import glob
from Bio import SeqIO

import os.path
import csv


folder = 'data/Legionella_pneumophila'
out_file = 'assemblies_chrs.csv'

ans = [['assembly', 'chromosome', 'description']]
unique_contigs = set()
unique_assembly_count = 0

for file in glob(f'{folder}/fna/*.fna'):
    if file.endswith('merged.fna'): continue
    unique_assembly_count += 1

    contigs = [contig for contig in SeqIO.parse(open(file), 'fasta') if not 'plasmid' in contig.description]
    assert len(contigs) == 1

    unique_contigs.add(contigs[0].id)

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
print("Unique contigs:", len(list(unique_contigs)))