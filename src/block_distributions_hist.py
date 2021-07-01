import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from src.utils.infercars_tools import parse_to_df, filter_unique_gene, dist_between_blocks, number_of_genomes_dist, \
    blocks_length_dist

sns.set(style="whitegrid", font="serif")
df = parse_to_df('data/E_coli/sibeliaz_out/fine/1000/blocks_coords.infercars')


def length_before():
    # df = filter_unique_gene(df)
    ds = dist_between_blocks(df)

    sns.distplot(ds, kde=False, bins=100, hist_kws={'log': True, 'alpha': 0.7})

    state = 'before'
    plt.ylabel('Number of blocks')
    plt.xlabel('Length in nucleotides, kilobases')
    # plt.title(f'Length between blocks distribution ({state} filtering)')

    plt.xlim(xmin=0)
    plt.tight_layout()
    plt.xticks(np.arange(0, 90000, 2e4), ('0', '20', '40', '60', '80'))
    plt.savefig(f'05_lengths_{state}_filtering.pdf')
    plt.show()


def length_after():
    df_filtered = filter_unique_gene(df)
    ds = dist_between_blocks(df_filtered)

    sns.distplot([])
    sns.distplot(ds, kde=False, bins=100, hist_kws={'log': True, 'alpha': 0.7})

    state = 'after'
    plt.ylabel('Number of blocks')
    plt.xlabel('Length in nucleotides, kilobases')
    # plt.title(f'Length between blocks distribution ({state} filtering)')

    plt.tight_layout()
    plt.xlim(xmin=0)
    plt.xticks(np.arange(0, 1.2e6, 2e5), ('0', '200', '400', '600', '800', '1000'))
    plt.savefig(f'06_lengths_{state}_filtering.pdf')
    plt.show()


def block_length():
    df_filtered = filter_unique_gene(df)
    ds_before = blocks_length_dist(df)
    ds_after = blocks_length_dist(df_filtered)

    bins = np.linspace(min(ds_before), max(ds_before), 100)
    sns.distplot(ds_before, kde=False, bins=bins, hist_kws={'log': True, 'alpha': 0.7}, label='not-common')
    sns.distplot(ds_after, kde=False, bins=bins, hist_kws={'log': True, 'alpha': 0.7}, label='common')

    state = 'after'
    plt.ylabel('Number of blocks')
    plt.xlabel('Length in nucleotides')
    plt.xlim(xmin=0)
    # plt.title(f'Blocks length distribution ({state} filtering)')

    plt.legend(loc=1)

    plt.tight_layout()
    plt.savefig(f'01_block_lengths_log.pdf')
    plt.show()


def block_length_example():
    df_filtered = filter_unique_gene(df)
    ds_before = blocks_length_dist(df)
    ds_after = blocks_length_dist(df_filtered)

    bins = np.linspace(min(ds_before), max(ds_before), 15)
    sns.distplot(ds_before, kde=False, bins=bins, hist_kws={'log': True}, label='not-common')
    sns.distplot(ds_after, kde=False, bins=bins, hist_kws={'log': True}, label='common')

    plt.ylabel('Number of blocks')
    plt.xlabel('Length in nucleotides')
    plt.xlim(xmin=0)
    plt.title(f'Blocks length distribution (same bins)')

    plt.legend(loc=1)

    plt.tight_layout()
    plt.savefig(f'same_bins.pdf')
    plt.show()


def number_of_genomes():
    vs = number_of_genomes_dist(df)
    sns.distplot(vs, kde=False, bins=100, hist_kws={'log': True, 'alpha': 0.7})

    plt.ylabel('Number of blocks')
    plt.xlabel('Number of genomes')
    # plt.title(f'Frequency of blocks within the 414 genomes')

    plt.xlim(xmin=0, xmax=414)
    plt.tight_layout()
    plt.savefig(f'03_blocks_frequency_hist.pdf')
    plt.show()


def number_of_genomes_weighted():
    vs = []
    ws = []
    for _, df_block in df.groupby('block'):
        vs.append(len(df_block.species.unique()))
        lens = [row['chr_end'] - row['chr_beg'] for _, row in df_block.iterrows()]
        ws.append(np.mean(lens))

    sns.distplot(vs, kde=False, bins=100, hist_kws={'log': True, 'alpha': 0.7, 'weights': ws})

    plt.ylabel('Length of fragments that are present\n in n genomes, nucleotides')
    plt.xlabel('Number of genomes')
    # plt.title(f'Frequency of blocks within the 414 genomes')

    plt.xlim(xmin=0, xmax=414)
    plt.tight_layout()
    plt.savefig(f'04_blocks_frequency_hist_weighted.pdf')
    plt.show()


def scatter_len_genomes_count():
    xs = []
    ys = []
    for block, df_block in df.groupby('block'):
        xs.append(len(df_block.species.unique()))
        lens = [row['chr_end'] - row['chr_beg'] for _, row in df_block.iterrows()]
        ys.append(np.mean(lens))
        print(block, len(df_block.species.unique()), round(np.mean(lens)), sep=',')

    plt.scatter(xs, ys, s=1)

    plt.xlabel('Number of genomes')
    plt.ylabel('Length of blocks')
    # plt.title(f'Dependency of number of genomes and length')

    plt.xlim(xmin=0 - 3, xmax=max(xs) + 3)
    plt.ylim(ymin=0)

    plt.subplots_adjust(left=0.14, right=0.99)

    plt.tight_layout()
    plt.savefig('02_scatter_number_length.pdf')
    plt.show()


def new_blocks(permutations=6543):
    block_sets = [set(df_sp.block.unique()) for _, df_sp in df.groupby('species')]

    oss = []
    for _ in range(permutations):
        block_sets = np.random.permutation(block_sets)
        os = []
        accumulate_set = set()
        for bs in block_sets:
            left = bs - accumulate_set
            os.append(len(left))
            accumulate_set |= bs

        oss.append(os)

    oss = np.array(oss)
    xs = list(range(1, len(block_sets) + 1))

    plt.plot(xs, np.median(oss, axis=0), label='median (different permutations)')
    plt.fill_between(xs, np.percentile(oss, 5, axis=0), np.percentile(oss, 95, axis=0), alpha=0.4,
                     label='90% confidence interval')

    plt.xlabel('number of genomes')
    plt.ylabel('new blocks')
    plt.legend(loc='upper right')

    plt.ylim(ymax=np.percentile(oss, 95, axis=0)[1], ymin=0)
    plt.xlim(xmin=0, xmax=414)

    # plt.subplots_adjust(top=0.99, right=0.99)
    plt.tight_layout()
    plt.savefig('08_new_blocks.pdf')
    plt.show()


def pan_blocks(permutations=3000):
    block_sets = [set(df_sp.block.unique()) for _, df_sp in df.groupby('species')]

    oss = []
    for _ in range(permutations):
        block_sets = np.random.permutation(block_sets)
        os = []
        accumulate_set = set()
        for bs in block_sets:
            accumulate_set |= bs
            os.append(len(accumulate_set))
        oss.append(os)

    oss = np.array(oss)
    xs = list(range(1, len(block_sets) + 1))

    plt.plot(xs, np.median(oss, axis=0), label='median (different permutations)', color='indianred')
    plt.fill_between(xs, np.percentile(oss, 5, axis=0), np.percentile(oss, 95, axis=0), alpha=0.4,
                     label='90% confidence interval', color='indianred')

    plt.xlabel('number of genomes')
    plt.ylabel('Pan-blocks count')
    plt.legend(loc='lower right')

    plt.ylim(ymin=np.percentile(oss, 5, axis=0)[0])
    plt.xlim(xmin=0, xmax=414)
    # plt.subplots_adjust(top=0.99, right=0.99)
    plt.tight_layout()
    plt.savefig('07_pan_blocks.pdf')
    plt.show()


length_before()
