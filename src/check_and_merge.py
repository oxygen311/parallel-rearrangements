import glob
from Bio import SeqIO

folder = 'data/Staphylococcus_aureus/fna'

all_contigs = []
for file in glob.glob(f'{folder}/*.fna'):
    if 'merged' in file: continue
    print(file)

    contigs = [contig for contig in SeqIO.parse(open(file), 'fasta') if not 'plasmid' in contig.description]
    assert len(contigs) == 1

    # if len(all_contigs) > 100:
    #     break
    all_contigs += contigs

# SeqIO.write(all_contigs, open(f'{folder}/merged_100.fna', 'w'), 'fasta')
SeqIO.write(all_contigs, open(f'{folder}/merged.fna', 'w'), 'fasta')