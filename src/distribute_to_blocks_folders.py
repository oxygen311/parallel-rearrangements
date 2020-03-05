import os

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pprint import pprint
from collections import defaultdict


folder = 'data/Bacillus_anthracis/'
sibeliaz_folder = folder + 'sibeliaz/fine/5000/'
blocks_faa_folder = folder + 'faa_from_blocks/'
csv_file = 'blocks_coords_unique_gene.csv'
assemblies_chr_file = 'assemblies_chrs.csv'

df = pd.read_csv(sibeliaz_folder + csv_file)


def make_chr_to_assembly_dict():
    df = pd.read_csv(folder + assemblies_chr_file)
    d = dict()
    for i, row in df.iterrows():
        d[row['chromosome']] = row['assembly']
    return d

chr_to_assembly = make_chr_to_assembly_dict()
used_ids_count = 0

for seq, df_seq in df.groupby('Seq'):
    print(seq)
    assembly = chr_to_assembly[seq]
    records = SeqIO.parse(folder + 'gbff/' + assembly + '.gbff', 'genbank')

    for record in records:
        # filtering plasmids
        if 'plasmid' in record.description: continue

        seqs = []
        for i, row_block in df_seq.iterrows():
            block_start = row_block['Start']
            block_end = row_block['End']

            if row_block['Strand'] == '+':
                sub_record = record[block_start:block_end]
            else:
                sub_record = record[block_end:block_start]

            used_ids = set()
            for feature in sub_record.features:
                if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                    q = feature.qualifiers
                    protein_id = q['protein_id'][0]
                    translation = q['translation'][0]

                    assert len(q['translation']) == 1
                    if protein_id in used_ids:
                        used_ids_count += 1
                        continue

                    simple_seq_r = SeqRecord(
                        Seq(translation),
                        'lcl|' + record.id + '_prot_' + protein_id,
                        description=f'[locus_tag={q["locus_tag"][0]}] [protein={q["product"][0]}] [protein_id={q["protein_id"][0]}] [location={feature.location}] [gbkey=CDS]')

                    used_ids.add(protein_id)
                    seqs.append(simple_seq_r)

        cur_folder = blocks_faa_folder
        if not os.path.exists(cur_folder): os.makedirs(cur_folder)
        SeqIO.write(seqs, cur_folder + f'{record.id}.faa', 'fasta')

print('Used ids count:', used_ids_count)