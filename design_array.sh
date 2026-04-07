#!/bin/bash

#SBATCH --job-name=amplicon
#SBATCH -p youlab-gpu
#SBATCH --time=24:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

#SBATCH --output=slurm_outputs/amplicon-%A_%a.out
#SBATCH --error=slurm_outputs/amplicon-%A_%a.err

#SBATCH --array=1-20        # number of random seeds; more seeds = better chance
                            # of finding the maximum-coverage solution

source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate amplicon_design
cd ~/xiaoli_PCR_design
export PYTHONPATH=~/xiaoli_PCR_design:$PYTHONPATH
mkdir -p slurm_outputs

# ── Edit these paths/lengths to match kmer_generation.py and candidates.sh ─
RANGES_CSV="./barcode_ranges_in_genome.csv"
ONLY_GENOME_DIR="./only_genome"
BACKGROUND_DIR="./genome_and_plasmids_within_host"
MIN_LENGTH=20
MAX_LENGTH=20      # set equal to MIN_LENGTH for a single length, or different for a range
# ──────────────────────────────────────────────────────────────────────────

SEED=${SLURM_ARRAY_TASK_ID}

echo "Running seed ${SEED} on $(hostname)"
python design_amplicon.py "$SEED" "$RANGES_CSV" "$ONLY_GENOME_DIR" \
                           "$BACKGROUND_DIR" "$MIN_LENGTH" "$MAX_LENGTH"
