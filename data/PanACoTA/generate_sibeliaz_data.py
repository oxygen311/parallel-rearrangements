import os

import pandas as pd

from Bio import SeqIO
from glob import glob

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--folder', '-f', required=True, help='Folder with PanACoTA output.')
args = parser.parse_args()
folder = args.folder if args.folder[-1] == '/' else args.folder + '/'

gembase_file = glob(folder + '2-annotate_module/LSTINFO-LSTINFO*')[0]
output_file = folder + '7-sibeliaz/for_sibeliaz.fna'
os.makedirs(folder + '7-sibeliaz/', exist_ok=True)

gembase_df = pd.read_csv(gembase_file, sep='\t')
all_contigs = []

for _, row in gembase_df.iterrows():
    print(row['gembase_name'], row['orig_name'])

    contigs = [contig for contig in SeqIO.parse(open(row['orig_name']), 'fasta') if not 'plasmid' in contig.description]

    if len(contigs) > 1:
        print('WARNING: Skip', row['gembase_name'])
        continue
    assert len(contigs) == 1

    contigs[0].id = row['gembase_name']

    all_contigs += contigs

SeqIO.write(all_contigs, open(output_file, 'w'), 'fasta')