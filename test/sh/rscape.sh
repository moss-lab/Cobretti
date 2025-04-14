#!/bin/bash -l
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=rscape
#SBATCH --mail-user=ecoppen@iastate.edu
#SBATCH --mail-type=END,FAIL

module load py-biopython
module load gcc
module load ghostscript
module load gnuplot

for f in *.stockholm; do /lustre/hdd/LAS/wmoss-lab/programs/rscape_v2.0.0.k/bin/R-scape -s --ntree 10 $f; done
sed -i.bak "/#=GF R2R*/d" *.sto
for g in *.sto; do /lustre/hdd/LAS/wmoss-lab/programs/R2R-1.0.6/src/r2r --disable-usage-warning $g $(basename $g sto)pdf; done
gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=All.Rscape.pdf -dBATCH *.R2R.pdf
wait;
python /lustre/hdd/LAS/wmoss-lab/ecoppen/benchmarks/cobretti_nova/Cobretti/cobretti.py -stage 1CA -email ecoppen@iastate.edu &
wait;
