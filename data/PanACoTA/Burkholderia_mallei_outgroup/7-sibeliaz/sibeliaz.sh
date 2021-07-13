#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=1-0
sibeliaz -k 15 -b 300 -t 32 -a 640 -n -o sibeliaz_out for_sibeliaz_2_contigs.fna
