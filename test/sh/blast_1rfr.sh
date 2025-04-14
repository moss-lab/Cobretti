#!/bin/bash -l
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=blast_1rfr
#SBATCH --mail-user=ecoppen@iastate.edu
#SBATCH --mail-type=END,FAIL

module load ncbi-rmblastn
echo ${NCBI_BLAST_DB_PATH}: /lustre/hdd/LAS/BioDatabase/ncbi/blast-db/latest
module load blast-plus
export BLASTDB=/lustre/hdd/LAS/BioDatabase/ncbi/blast-db/latest

blastn -task blastn -db nt -query /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1rfr.fasta -out /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1rfr_db_raw.txt -outfmt "6 sacc sseq" -max_target_seqs 2500 -num_threads 8 -max_hsps 1;
wait;