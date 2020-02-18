import os

import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pprint import pprint
from collections import defaultdict


folder = 'data/Bacillus_anthracis/'
sibeliaz_folder = folder + 'sibeliaz/fine/5000/'
blocks_faa_folder = folder + 'faa_blocks/'
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

for seq, df_seq in df.groupby('Seq'):
    print(seq)
    assembly = chr_to_assembly[seq]
    records = SeqIO.parse(folder + 'gbff/' + assembly + '.gbff', 'genbank')

    for record in records:
        # filtering plasmids
        if 'plasmid' in record.description: continue

        for i, row_block in df_seq.iterrows():
            seqs = []
            block = row_block['Block']
            block_start = row_block['Start']
            block_end = row_block['End']

            if row_block['Strand'] == '+':
                sub_record = record[block_start:block_end]
            else:
                sub_record = record[block_end:block_start]

            for feature in sub_record.features:
                if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                    q = feature.qualifiers

                    assert len(q['translation']) == 1
                    simple_seq_r = SeqRecord(
                        Seq(q['translation'][0]),
                        'lcl|' + record.id + '_prot_' + q['protein_id'][0],
                        description=f'[locus_tag={q["locus_tag"][0]}] [protein={q["product"][0]}] [protein_id={q["protein_id"][0]}] [location={feature.location}] [gbkey=CDS]')

                    seqs.append(simple_seq_r)

            cur_folder = blocks_faa_folder + f'block_{block}/'
            if not os.path.exists(cur_folder): os.makedirs(cur_folder)
            SeqIO.write(seqs, cur_folder + f'{record.id}.faa', 'fasta')