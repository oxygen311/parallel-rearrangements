from collections import defaultdict

import numpy as np

def get_class_colors(labels_dict, genome_color):
    genome_classes_colors = defaultdict(list)
    for genome, color in genome_color.items():
        splitted = labels_dict[genome].split()
        genome_class = splitted[0] + " " + splitted[1]

        genome_classes_colors[genome_class].append(color)
    return genome_classes_colors

def white_proportion(colors):
    return np.mean(list(map(lambda c: c == 0, colors)))

def count_shigella_differs_from_white(class_colors):
    return sum(
        white_proportion(class_colors[cls]) < 1
        for cls in
        class_colors.keys() - {'Escherichia coli'})

def count_shigella_differs_from_value(class_colors, value, threshhold=0.3):
    sh_keys = class_colors.keys() - {'Escherichia coli'}
    return sum(np.mean(class_colors[cls]) - value >= threshhold for cls in sh_keys), \
           sum(np.mean(class_colors[cls]) - value < threshhold for cls in sh_keys)