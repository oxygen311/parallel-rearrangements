import os

import pandas as pd

from Bio import SeqIO
from glob import glob

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--folder', '-f', required=True, help='Folder with PanACoTA output.')
parser.add_argument('--contigs', '-c', required=True, help='Number of maximum contigs to take from every genome.',
                    type=int)
args = parser.parse_args()

folder = args.folder if args.folder[-1] == '/' else args.folder + '/'

gembase_file = glob(folder + '2-annotate_module/LSTINFO-LSTINFO*')[0]
output_file = folder + f'7-sibeliaz/for_sibeliaz_{args.contigs}_contigs.fna'
os.makedirs(folder + '7-sibeliaz/', exist_ok=True)

gembase_df = pd.read_csv(gembase_file, sep='\t')
all_contigs = []
from collections import Counter
cnt = Counter()

for _, row in gembase_df.iterrows():
    print(row['gembase_name'], row['orig_name'])

    contigs = [contig for contig in SeqIO.parse(open(row['orig_name']), 'fasta')][0:args.contigs]

    for i, contig in enumerate(contigs):
        # print(contig.id, contig.name, contig.description)
        contigs[i].id = row['gembase_name'] + '.' + str(i + 1) + ' ' + contig.description

    print(len(contigs))
    # assert len(contigs) == args.contigs
    cnt[len(contigs)] += 1

    all_contigs += contigs

print(cnt)
SeqIO.write(all_contigs, open(output_file, 'w'), 'fasta')