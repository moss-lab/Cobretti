#!/bin/bash -l
#SBATCH --nodelist=whatever
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=cobretti
#SBATCH --mail-user=[example]@iastate.edu
#SBATCH --mail-type=ALL

module load py-biopython/1.70-py3-wos466g
module load python/3.6.5-fwk5uaj

python cobretti.py -stage 1A -email [example]@iastate.edu
