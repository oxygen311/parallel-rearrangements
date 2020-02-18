import glob
from Bio import SeqIO

folder = 'data/Legionella_pneumophila/fna'

all_contigs = []
for file in glob.glob(f'{folder}/*.fna'):
    if file.endswith('merged.fna'): continue
    print(file)

    contigs = [contig for contig in SeqIO.parse(open(file), 'fasta') if not 'plasmid' in contig.description]
    assert len(contigs) == 1

    all_contigs += contigs

SeqIO.write(all_contigs, open(f'{folder}/merged.fna', 'w'), 'fasta')