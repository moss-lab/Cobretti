#!/bin/bash -l
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=scanfold2_1a51
#SBATCH --mail-user=ecoppen@iastate.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --export=NONE

module purge
module load micromamba

micromamba activate /lustre/hdd/LAS/wmoss-lab/programs/envs/ScanFold2
wait;
python /lustre/hdd/LAS/wmoss-lab/scripts/ScanFold2.0-inforna/ScanFoldBothForInforna.py /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta --folder_name 1a51 --global_refold &
wait;
