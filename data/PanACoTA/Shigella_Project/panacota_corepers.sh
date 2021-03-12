#!/bin/bash -i
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-0
#SBATCH -o panacota_corepers.out
#SBATCH -e panacota_corepers.err
conda activate panacota
PanACoTA corepers -p 3-pangenome_module/PanGenome-ESCO414.All.prt-clust-0.95-mode1-th32_2021-03-12_22-33-08.tsv.lst -o 4-corepers_module_0.95
