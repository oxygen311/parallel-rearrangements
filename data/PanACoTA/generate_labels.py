import pandas as pd

from Bio import SeqIO
from glob import glob

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--folder', '-f', required=True, help='Folder with PanACoTA output.')
args = parser.parse_args()
folder = args.folder if args.folder[-1] == '/' else args.folder + '/'

summary_file = glob(folder + '1-prepare_module/assembly_summary*')[0]
summary_df = pd.read_csv(summary_file, sep='\t')

gembase_file = glob(folder + '2-annotate_module/LSTINFO-LSTINFO*')[0]
gembase_df = pd.read_csv(gembase_file, sep='\t')

outfile = folder + 'labels.csv'

# print(summary_df.columns)

gembase_df['file'] = gembase_df['orig_name'].str.split('/').str[-1]

gembase_df['assembly_accession'] = gembase_df['file'].str.split('_').str[0] + '_' + \
                                   gembase_df['file'].str.split('_').str[1]
gembase_df['asm_name'] = gembase_df['file'].str.split('_', 2).str[2].str.rsplit('_', 2).str[0]

# merged = pd.merge(summary_df, gembase_df, on=['assembly_accession', 'asm_name'])
merged = pd.merge(summary_df, gembase_df, on=['assembly_accession'])

assert len(merged) == len(gembase_df)

# print(merged['infraspecific_name'].str.replace('strain=', ''))
# print(merged["infraspecific_name"].str.contains(merged["infraspecific_name"]))

merged['name_no_strain'] = merged['infraspecific_name'].str.replace('strain=', '').astype(str)
merged['label'] = merged.apply(
    lambda x: x.organism_name if x.name_no_strain in x.organism_name else x.organism_name + ' ' + x.name_no_strain,
    axis=1)

merged.rename(columns={'gembase_name': 'strain'}, inplace=True)
merged.to_csv(outfile, index=False)