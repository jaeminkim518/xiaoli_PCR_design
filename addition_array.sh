#!/bin/bash

#SBATCH --job-name=add_amp
#SBATCH -p youlab-gpu
#SBATCH --time=24:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

#SBATCH --output=slurm_outputs/addition-%A_%a.out
#SBATCH --error=slurm_outputs/addition-%A_%a.err

#SBATCH --array=1-20        # number of random seeds; more seeds = better chance
                            # of placing every new target

source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate amplicon_design
cd ~/xiaoli_PCR_design
export PYTHONPATH=~/xiaoli_PCR_design:$PYTHONPATH
mkdir -p slurm_outputs

# ── Edit these to match kmer_generation.py, candidates.sh and design_array.sh ─
ONLY_GENOME_DIR="./only_genome"
BACKGROUND_DIR="./genome_and_plasmids_within_host"
MIN_LENGTH=20
MAX_LENGTH=20      # set equal to MIN_LENGTH for a single length

# The panel to keep fixed — a validated_primers CSV from Stage 3
LOCKED_CSV="./output/barcodes_k20_amplicon_design/validated_primers_k20-seed7.csv"

# Optional. Leave empty to design every isolate/region not already locked.
TARGETS=""         # file of 'SEQ_ID,Region' or bare SEQ_ID lines
UNLOCK=""          # file of locked targets to release and re-design
EXTRA_ARGS=""      # e.g. "--retry" or "--skip-locked-pair-audit"
# ──────────────────────────────────────────────────────────────────────────

SEED=${SLURM_ARRAY_TASK_ID}

# The audit does not depend on the seed, so run it once (task 1) and skip it
# on the rest to save wall time.
if [[ "${SLURM_ARRAY_TASK_ID}" -ne 1 ]]; then
    EXTRA_ARGS="${EXTRA_ARGS} --skip-locked-audit"
fi

ARGS=(--locked "$LOCKED_CSV"
      --only-genome "$ONLY_GENOME_DIR"
      --background "$BACKGROUND_DIR"
      --length "$MIN_LENGTH" --max-length "$MAX_LENGTH"
      --seed "$SEED")

[[ -n "$TARGETS" ]] && ARGS+=(--targets "$TARGETS")
[[ -n "$UNLOCK"  ]] && ARGS+=(--unlock  "$UNLOCK")
# shellcheck disable=SC2206
[[ -n "$EXTRA_ARGS" ]] && ARGS+=($EXTRA_ARGS)

echo "Running seed ${SEED} on $(hostname)"
echo "Locked panel: ${LOCKED_CSV}"

python add_amplicons.py "${ARGS[@]}"
