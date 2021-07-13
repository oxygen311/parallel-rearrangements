#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=1-0
maf2synteny -s fine.txt -o fine -b 5000 blocks_coords.gff
maf2synteny -s fine.txt -o fine -b 3000 blocks_coords.gff
maf2synteny -s fine_500.txt -o fine -b 1000 blocks_coords.gff
