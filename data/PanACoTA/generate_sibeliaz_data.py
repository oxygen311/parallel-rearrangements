import pandas as pd

from Bio import SeqIO

gembase_file = 'Streptococcus_pyogenes/2-annotate_module/LSTINFO-LSTINFO-1314-filtered-1e-06_0.06.lst'
output_file = 'Streptococcus_pyogenes/7-sibeliaz/for_sibeliaz.fna'

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