#!/bin/bash -l
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=cmbuilder_1
#SBATCH --mail-user=ecoppen@iastate.edu
#SBATCH --mail-type=FAIL

module load perl
module load infernal

export PERL5LIB=/lustre/hdd/LAS/wmoss-lab/programs/lib/perl5
export PERL5LIB=/lustre/hdd/LAS/wmoss-lab/programs/RNAFramework/lib
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_2_hfold_scanfold.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_1_hfold.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_3_hfold.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_2_hfold.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_1_knotty_2.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_2_knotty_2.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_3_knotty_2.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_1_knotty_1.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_2_knotty_1.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
perl /lustre/hdd/LAS/wmoss-lab/scripts/labtools/cm-builder -s /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/sequences/1a51.fasta -m /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/pk_motifs/1a51_motif_3_knotty_1.dbn -d /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/test/databases/1a51_db.fa -c 4 -T ./tmp1 -k -t 1 &
wait;
