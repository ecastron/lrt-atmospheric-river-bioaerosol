#!/bin/bash
#SBATCH --job-name=coxL_iqtree_01OCT26
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=05_tree_%j.log

# 05_tree.sh
# Maximum-likelihood tree of the motif-labeled CoxL references, the outgroup and the 89 dereplicated air candidates
# (MAFFT --auto alignment; IQ-TREE LG+G4 with 1,000 ultrafast bootstraps).
# Set COXL_DIR to this folder (default: current directory).
# Run    : sbatch 05_tree.sh

set -euo pipefail
cd "${COXL_DIR:-$PWD}"
set +u; module purge; module load MAFFT/MAFFT; set -u
mafft --version
mafft --auto --thread "$SLURM_CPUS_PER_TASK" tree_input.faa > tree_input.aln 2>/dev/null
set +u; module purge; module load Iqtree/iqtree; set -u
IQ=$(command -v iqtree2 || command -v iqtree || command -v iqtree3)
"$IQ" --version | head -1
"$IQ" -s tree_input.aln -m LG+G4 -B 1000 -T "$SLURM_CPUS_PER_TASK" --prefix coxL_iqtree -redo -quiet
echo done
