#!/bin/bash

#SBATCH -J cand_amp
#SBATCH -p youlab-gpu
#SBATCH -t 06:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

#SBATCH --array=0-17        # 0 to (number of isolates - 1); check with:
                            #   ls ./only_genome/*.fasta | grep -v '\._' | wc -l
#SBATCH -o slurm_outputs/%x_%A_%a.out
#SBATCH -e slurm_outputs/%x_%A_%a.err

source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate amplicon_design
cd /hpc/home/jk630/xiaoli_PCR_design
export PYTHONPATH="/hpc/home/jk630/xiaoli_PCR_design:${PYTHONPATH:-}"
mkdir -p slurm_outputs

# ── Edit these paths/lengths to match your setup ──────────────────────────
RANGES_CSV="./barcode_ranges_in_genome.csv"
ONLY_GENOME_DIR="./only_genome"
BACKGROUND_DIR="./genome_and_plasmids_within_host"
MIN_LENGTH=20
MAX_LENGTH=20      # set equal to MIN_LENGTH for a single length, or different for a range
# ──────────────────────────────────────────────────────────────────────────

SCRIPT=individual_target_amplicon_candidates.py

mapfile -d '' -t FASTAS < <(
    find "$ONLY_GENOME_DIR" -maxdepth 1 -type f -name '*.fasta' \
        ! -name '._*' -print0 | sort -z
)
N=${#FASTAS[@]}

if [[ "$SLURM_ARRAY_TASK_ID" -ge "$N" ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} >= N=${N}; exiting."
    exit 0
fi

FASTA="${FASTAS[$SLURM_ARRAY_TASK_ID]}"
SEQ_ID=$(basename "$FASTA" .fasta)

echo "Running SEQ_ID=$SEQ_ID on $(hostname)"
python "$SCRIPT" "$RANGES_CSV" "$ONLY_GENOME_DIR" "$BACKGROUND_DIR" \
                 "$SEQ_ID" "$MIN_LENGTH" "$MAX_LENGTH"
