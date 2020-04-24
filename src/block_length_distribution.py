from src.utils.infercars_tools import parse_to_df, filter_unique_gene, dist_between_blocks

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid", font="serif")

df = parse_to_df('data/E_coli/sibeliaz_out/fine/1000/blocks_coords.infercars')
# df = filter_unique_gene(df)

ds = dist_between_blocks(df)

sns.distplot(ds, kde=False, bins=100, hist_kws={'log': True})

state = 'before'
plt.ylabel('count')
plt.xlabel('length in nucleotides')
plt.title(f'Length between blocks distribution ({state} filtering)')

plt.savefig(f'lengths_{state}_filtering.pdf')
plt.show()

