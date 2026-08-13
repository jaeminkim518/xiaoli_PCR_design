#!/bin/bash
#SBATCH --job-name=amp_lockset
#SBATCH -p youlab-gpu
#SBATCH --time=24:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

#SBATCH --output=slurm_outputs/lockset-%A_%a.out
#SBATCH --error=slurm_outputs/lockset-%A_%a.err

#SBATCH --array=1-24

# ---------------------------------------------------------------------------
# Works for both --mode extend and --mode rescue.
#
# Submit with:
#   sbatch --export=ALL,MODE=extend,TARGET_COMM=comm3,KMER_SIZE=20,MIN_LEN=80,MAX_LEN=120,\
# LOCKED_CSV=./output/comm3_k20_L80-120_amplicon_design/validated_primers_k20-seed1.csv,\
# NEW_IDS=./new_isolates.txt rescue_array.sh
#
# Every variable below falls back to a default if it is not exported, so
# `set -u` cannot kill the job the way it does when a bare ${KMER_SIZE} is
# referenced but never passed.
# ---------------------------------------------------------------------------

source /hpc/group/youlab/zz294/miniconda3/etc/profile.d/conda.sh
conda activate amplicon_design

set -euo pipefail
mkdir -p slurm_outputs

: "${MODE:=extend}"
: "${TARGET_COMM:=comm3}"
: "${KMER_SIZE:=20}"
: "${MIN_LEN:=80}"
: "${MAX_LEN:=120}"
: "${LOCKED_CSV:?LOCKED_CSV must be exported (path to validated_primers CSV)}"
: "${NEW_IDS:=}"       # optional: file of newly added seq_ids (scopes the audit)
: "${TARGETS:=}"       # required for MODE=rescue: file of seq_ids to design
: "${EXTRA_ARGS:=}"    # e.g. "--retry" or "--skip-locked-pair-audit"

SEED=${SLURM_ARRAY_TASK_ID}

# The locked-set audit is identical for every seed. Run it only on task 1;
# the rest skip it to save wall time.
if [[ "${SLURM_ARRAY_TASK_ID}" -ne 1 ]]; then
    EXTRA_ARGS="${EXTRA_ARGS} --skip-locked-audit"
fi

ARGS=(--mode "${MODE}" --comm "${TARGET_COMM}" --kmer "${KMER_SIZE}"
      --min-len "${MIN_LEN}" --max-len "${MAX_LEN}"
      --seed "${SEED}" --locked "${LOCKED_CSV}")

[[ -n "${NEW_IDS}" ]] && ARGS+=(--new-ids "${NEW_IDS}")
[[ -n "${TARGETS}" ]] && ARGS+=(--targets "${TARGETS}")
# shellcheck disable=SC2206
[[ -n "${EXTRA_ARGS}" ]] && ARGS+=(${EXTRA_ARGS})

echo "Task ${SLURM_ARRAY_TASK_ID} | mode=${MODE} seed=${SEED} comm=${TARGET_COMM} k=${KMER_SIZE} L=${MIN_LEN}-${MAX_LEN}"
echo "Locked: ${LOCKED_CSV}"

python rescue_amplicon.py "${ARGS[@]}"
