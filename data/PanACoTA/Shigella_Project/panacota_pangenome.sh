#!/bin/bash -i
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-0
#SBATCH -o panacota_pangenome.out
#SBATCH -e panacota_pangenome.err
#SBATCH --nodelist=orthrus-[2]
conda activate panacota
PanACoTA pangenome -l 2-annotate_module/LSTINFO-list.lst -n ESCO414 -d 2-annotate_module/Proteins/ -o 3-pangenome_module -i 0.95 --threads 32
