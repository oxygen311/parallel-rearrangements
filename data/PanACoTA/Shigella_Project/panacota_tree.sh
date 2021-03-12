#!/bin/bash -i
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-0
#SBATCH -o panacota_tree.out
#SBATCH -e panacota_tree.err
conda activate panacota
PanACoTA tree -a 5-align_module_0.95/Phylo-ESCO414-0.95/ESCO414-0.95.grp.aln -o 6-tree_module_0.95 --boot 1000 --threads 32
