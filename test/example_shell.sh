#!/bin/bash -l
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=cobretti
#SBATCH --mail-user=ecoppen@iastate.edu
#SBATCH --mail-type=ALL

module load py-biopython
module load python

python ../cobretti.py -stage 1B -email ecoppen@iastate.edu
