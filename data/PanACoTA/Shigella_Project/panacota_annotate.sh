#!/bin/bash -i
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-0
#SBATCH -o panacota_annotate.out
#SBATCH -e panacota_annotate.err
conda activate panacota
PanACoTA annotate -d 1-chromosomes -r 2-annotate_module -n ESCO -l Chromosomes_fna/list.txt --cutn 0 --threads 32
