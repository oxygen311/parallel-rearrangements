from bg.tree import BGTree
from bg.multicolor import Multicolor
from bg.genome import BGGenome

from src.bp_analyze.tree_holder import TreeHolder
from src.bp_analyze.common import print_species_stats, make_labels_dict, get_genomes_contain_blocks, get_blocks_order
from src.bp_analyze.consistency_checker import TreeConsistencyChecker
from src.bp_analyze.shigella_common import get_class_colors, count_shigella_differs_from_value
from src.utils.infercars_tools import parse_to_df

from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import MDS
from collections import defaultdict

from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace


import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

tree_file = 'data/Streptococcus_pneumoniae/tree/RAxML_bestTree.concat_alignment.tree'

if __name__ == "__main__":
    tree = Tree(tree_file)
    prev_n = None
    for node in tree.traverse():
        if prev_n is not None:
            print(node.get_distance(prev_n))
        prev_n = node
    print(tree.get_distance("CP000920", "CM001835"))
    # print(tree)
