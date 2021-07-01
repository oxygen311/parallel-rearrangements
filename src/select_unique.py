import glob
from Bio import SeqIO
from bg.grimm import GRIMMReader
from collections import defaultdict

file = 'data/Staphylococcus_aureus/sibelia_out_a_1100_k_15_m_30/fine/1000/genomes_permutations.txt'
folder = 'data/Staphylococcus_aureus/fna'

def get_unique_contigs(file):
    with open(file) as f:
        l = f.readlines()

    seqs = defaultdict(list)

    for i in range(0, len(l), 2):
        name = GRIMMReader.parse_genome_declaration_string(l[i]).name
        data = GRIMMReader.parse_data_string(l[i + 1])[1]

        only_numbers = [d[1] for d in data]
        ind = only_numbers.index('1')
        data = data[ind:] + data[:ind]

        data_str = ' '.join(map(lambda a: f'{a[0]}{a[1]}', data))
        if data_str[0] == '-':
            data_str = data_str.replace('+', '?').replace('-', '+').replace('?', '-')

        seqs[data_str].append(name)

    return set(vs[0] for vs in seqs.values())

def select_contigs_fasta(folder, unique_contigs_ids):
    unique_contigs = []
    for file in glob.glob(f'{folder}/*.fna'):
        if 'merged' in file: continue
        unique_contigs.extend(contig for contig in SeqIO.parse(open(file), 'fasta') if contig.id in unique_contigs_ids)

    SeqIO.write(unique_contigs, open(f'{folder}/unique_merged.fna', 'w'), 'fasta')
    print(f'Wrote {len(unique_contigs)} unique contigs')

unique_contigs_ids = get_unique_contigs(file)
print(f'Got {len(list(unique_contigs_ids))} unique contigs')

select_contigs_fasta(folder, unique_contigs_ids)
