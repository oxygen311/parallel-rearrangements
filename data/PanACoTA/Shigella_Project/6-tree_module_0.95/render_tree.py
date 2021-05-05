from parebrick.tree.tree_holder import TreeHolder
from parebrick.utils.data.parsers import make_labels_dict

import logging

from collections import defaultdict


t = TreeHolder('asm_tree_reroot.tree', logging.getLogger(), labels_dict=make_labels_dict('LABELS.csv'))

t.count_innovations_fitch(defaultdict(int))

t.draw('test_2.svg', colors=['White'])