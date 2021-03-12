#!/bin/bash -i
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-0
#SBATCH -o panacota_align.out
#SBATCH -e panacota_align.err
conda activate panacota
PanACoTA align -c 4-corepers_module_0.95/PersGenome_PanGenome-ESCO414.All.prt-clust-0.95-mode1-th32_2021-03-12_22-33-08.tsv.lst_1.lst -l 2-annotate_module/LSTINFO-list.lst -n ESCO414-0.95 -d 2-annotate_module/ -o 5-align_module_0.95 --threads 32
