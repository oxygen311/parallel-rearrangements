from src.utils.infercars_tools import parse_to_df, export_df, filter_unique_gene

in_file = 'data/E_coli/sibeliaz_out/fine/1000/blocks_coords.infercars'
out_file = 'data/E_coli/sibeliaz_out/fine/1000/blocks_coords_unique_gene.infercars'

df = parse_to_df(in_file)
unique_df = filter_unique_gene(df)
export_df(unique_df, out_file)