from Bio import SeqIO
from glob import glob

import os.path

folder = 'data/EPEC/gbk/'
output_file = 'data/EPEC/merged.fna'

seqs = []
for file in glob(folder + '*.gb'):
    basename = os.path.basename(file)
    epecname = basename.split('.')[0]
    print(epecname)

    records = SeqIO.parse(file, "genbank")

    for i, seq_record in enumerate(records):
        seq_record.id = epecname + '.contig_' + str(i)
        seq_record.description = ''
        seqs.append(seq_record)


SeqIO.write(seqs, output_file, 'fasta')