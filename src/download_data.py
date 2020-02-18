from Bio import Entrez
from pprint import pprint
from collections import Counter

import pandas as pd

import os
import urllib.request
import urllib.error

Entrez.email = 'a.zabelkin@itmo.ru'
type = 'gbff'
term = 'Legionella pneumophila'

def get_assembly_summary(id):
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    return Entrez.read(esummary_handle)

def get_assemblies(term, download=True, base_folder='data/'):
    folder = base_folder + '_'.join(term.split()) + '/' + type
    if not os.path.exists(folder): os.makedirs(folder)

    handle = Entrez.esearch(db="assembly", term=term, retmax='100500')
    record = Entrez.read(handle)

    ids = record['IdList']
    print(f'found {len(ids)} ids')
    summaries = get_assembly_summary(','.join(ids))

    i = 0
    for summary in filter(lambda s: s['AssemblyStatus'] == 'Complete Genome',
                          summaries['DocumentSummarySet']['DocumentSummary']):
        print('Current:', summary['AssemblyName'])

        for sp in summary['Biosource']['InfraspeciesList']:
            # print(sp)
            print(sp['Sub_type'])
            assert sp['Sub_type'] == 'strain'
        # get ftp link
        url = summary['FtpPath_RefSeq']
        if url == '': continue
        label = os.path.basename(url)

        # get the fasta link - change this to get other formats
        link = os.path.join(url, label + f'_genomic.{type}.gz')
        file_path = f'{folder}/{label}.{type}.gz'
        if download and not os.path.exists(file_path):
            print('Downloading to:', file_path)
            urllib.request.urlretrieve(link, file_path)

        i += 1

    print(f'Got {i} genomes')

get_assemblies(term, True)
