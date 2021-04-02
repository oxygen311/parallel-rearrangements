import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import math
import numpy as np

origin_file = 'Streptococcus_pyogenes/ori_ters.csv'
lengths_file = 'Streptococcus_pyogenes/8-parebrick/3000/preprocessed_data/genomes_lengths.csv'
blocks_file = 'Streptococcus_pyogenes/8-parebrick/3000/preprocessed_data/blocks_coords.csv'
outfile = 'Streptococcus_pyogenes/distances_from_middle_of_breaks_to_ori.csv'

origin_df = pd.read_csv(origin_file, sep=',')
lengths_df = pd.read_csv(lengths_file, sep=',')
blocks_df = pd.read_csv(blocks_file, sep=',')

strain_to_ori = {row['gembase_name']: row['origin'] for _, row in origin_df.iterrows()}
strain_to_len = {row['Genome']: row['Length'] for _, row in lengths_df.iterrows()}

strain_to_breakpoints = {}
ds = []

for strain, df_strain in blocks_df.groupby('species'):
    len = strain_to_len[strain]
    ori = strain_to_ori[strain]

    bl1 = df_strain[df_strain['block'] == 22].iloc[0]
    bl2 = df_strain[df_strain['block'] == 124].iloc[0]

    br1 = bl1.chr_beg if bl1.orientation == '+' else bl1.chr_end
    br2 = bl2.chr_beg if bl2.orientation == '-' else bl2.chr_end

    strain_to_breakpoints[strain] = (br1, br2)

    middle1 = (br1 + br2) // 2
    middle2 = ((br1 + br2 + len) // 2) % len

    # print('middles:', middle1, middle2)

    dist1 = (ori - middle1 + len) % len
    dist2 = (middle1 - ori + len) % len
    dist3 = (ori - middle2 + len) % len
    dist4 = (middle2 - ori + len) % len

    min_dist = min(dist1, dist2, dist3, dist4)
    ds.append([strain, min_dist])

ds_df = pd.DataFrame(ds, columns=['strain', 'dist from middle to origin'])
ds_df.sort_values('dist from middle to origin', inplace=True)

print(ds_df)
ds_df.to_csv(outfile, index=False)
# sns.set_style('whitegrid')
# plt.hist(ds, bins=100)
# plt.show()