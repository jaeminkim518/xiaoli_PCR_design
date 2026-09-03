"""
Select an optimal combination of amplicons and primers across all isolates.

Usage:
    python design_amplicon.py <seed> <RANGES_CSV> <ONLY_GENOME_DIR> \\
                               <BACKGROUND_DIR> <BARCODE_LENGTH>
"""

import os
import sys
import random
# Ensure the script's own directory is on the path so local modules are found
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sys
import time as timer
import math

import numpy as np
import primer3

from Bio.Seq import Seq
from read_file_func import (
    parse_all_records_from_dir,
    load_all_candidate_primers,
)
from sequence_alignment import is_primer_pair_specific


# ---------------------------------------------------------------------------
# Hamming distance
# ---------------------------------------------------------------------------

def calculate_hamming_distance(s1, s2):
    """
    Returns Hamming distance for equal-length strings, or None if lengths differ.
    Callers should treat None as 'no constraint' (different lengths are always
    distinguishable regardless of sequence similarity).
    """
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


# ---------------------------------------------------------------------------
# Primer-dimer check
# ---------------------------------------------------------------------------

def has_cross_dimer(fwd1, rev1, fwd2, rev2,
                    temp_c=60.0, dv_conc=1.5,
                    max_dimer_dg=-9000.0, max_3prime_dg=-9000.0):
    for f, r in [(fwd1, rev2), (fwd2, rev1)]:
        if primer3.bindings.calc_heterodimer(f, r, temp_c=temp_c, dv_conc=dv_conc).dg < max_dimer_dg:
            return True
        if primer3.bindings.calc_end_stability(seq1=f, seq2=r,
                                                mv_conc=50.0, dv_conc=dv_conc,
                                                temp_c=temp_c).dg < max_3prime_dg:
            return True
    return False


# ---------------------------------------------------------------------------
# Weighted combination selection
# ---------------------------------------------------------------------------

def select_optimal_combination_weighted(all_candidates, min_hamming_distance):
    """
    Returns (barcode_counts, final_selection, failed_ids, failed_causes).
    Weighted ordering: isolates with fewer candidates are prioritised so they
    get the best chance to claim a barcode before the pool shrinks.
    """
    if not all_candidates:
        return {}, {}, list(all_candidates.keys()), []

    seq_ids       = list(all_candidates.keys())
    barcode_counts = {pid: len(all_candidates.get(pid, [])) for pid in seq_ids}
    max_sqrt      = math.sqrt(max(barcode_counts.values()) + 1) if barcode_counts else 1
    weights       = [max_sqrt - math.log(barcode_counts[pid] + 1) + 1 for pid in seq_ids]
    total_weight  = sum(weights)
    norm_weights  = [w / total_weight for w in weights]

    ordered_seqs = [seq_ids[i] for i in np.random.choice(
        len(seq_ids), size=len(seq_ids), replace=False, p=norm_weights
    )]

    final_selection = {}
    failed_ids      = []
    failed_causes   = []

    for seq_id in ordered_seqs:
        candidates_for_seq = list(all_candidates[seq_id])
        random.shuffle(candidates_for_seq)
        target_compatible = False
        hamming_issue = primer_dimer = 0

        for candidate in candidates_for_seq:
            cand_barcode = candidate['target_barcode']
            cand_fwd     = candidate['fwd_primer']
            cand_rev     = candidate['rev_primer']

            primer_ok = True
            for _, existing in final_selection.items():
                # Test 1: barcode Hamming distance (only for equal-length barcodes;
                # different lengths are inherently distinguishable)
                dist = calculate_hamming_distance(cand_barcode, existing['target_barcode'])
                if dist is not None and dist < min_hamming_distance:
                    primer_ok = False
                    hamming_issue += 1
                    break
                # Test 2: cross-dimer
                if has_cross_dimer(cand_fwd, cand_rev,
                                   existing['fwd_primer'], existing['rev_primer']):
                    primer_ok = False
                    primer_dimer += 1
                    break

            if primer_ok:
                target_compatible    = True
                final_selection[seq_id] = candidate
                break

        if not target_compatible:
            failed_ids.append(seq_id)
            failed_causes.append([hamming_issue, primer_dimer])

    assert set(final_selection) | set(failed_ids) == set(all_candidates), \
        "Some seq_ids are missing from both final_selection and failed_ids."
    return barcode_counts, final_selection, failed_ids, failed_causes


# ---------------------------------------------------------------------------
# Cross-reactivity check
# ---------------------------------------------------------------------------

def check_cross_reactivity(candidate_primers, all_sequences,
                            max_pcr_product=300, seed_len=4, max_mismatches=3):
    """
    All-against-all specificity check.
    all_sequences should include ALL genome + plasmid records for all isolates.
    Returns set of (seqA, seqB) pairs that cross-amplify.
    """
    problematic_pairs = set()
    primer_list = list(candidate_primers.items())

    for i in range(len(primer_list)):
        for j in range(i + 1, len(primer_list)):
            seq_A, primer_A = primer_list[i]
            seq_B, primer_B = primer_list[j]
            cross_pairs = [
                (primer_A['fwd_primer'], primer_B['fwd_primer']),
                (primer_A['fwd_primer'], primer_B['rev_primer']),
                (primer_A['rev_primer'], primer_B['fwd_primer']),
                (primer_A['rev_primer'], primer_B['rev_primer']),
            ]
            for fwd, rev in cross_pairs:
                if not is_primer_pair_specific(
                    fwd, rev, "",   # no target_id — checking cross-reactions
                    all_sequences,
                    max_pcr_product=max_pcr_product,
                    seed_len=seed_len,
                    max_mismatches=max_mismatches,
                ):
                    problematic_pairs.add((seq_A, seq_B))
                    break

    return problematic_pairs


def resolve_conflicts_greedy(problematic_pairs, candidate_counts):
    """Greedily remove highest-conflict nodes to maximise the retained set."""
    if not problematic_pairs:
        return set()

    graph: dict[str, set] = {}
    for a, b in problematic_pairs:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)

    alive     = set(graph.keys())
    failed_ids = set()

    def degree(n):
        return sum(1 for nbr in graph.get(n, set()) if nbr in alive)

    while True:
        candidates = [n for n in alive if degree(n) > 0]
        if not candidates:
            break
        pick = max(candidates, key=lambda n: (degree(n), -candidate_counts[n]))
        alive.remove(pick)
        failed_ids.add(pick)

    return failed_ids


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_primer_results(output_dir, final_primers, filename):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    with open(output_path, 'w', newline='') as f:
        header = ["Seq_ID", "Target_Barcode", "Forward_Primer", "Reverse_Primer",
                  "Product_Size", "Fwd_Tm", "Rev_Tm", "Amplicon_Sequence"]
        f.write(",".join(header) + "\n")
        for seq_id, data in sorted(final_primers.items()):
            row = [str(x) for x in [
                seq_id, data['target_barcode'], data['fwd_primer'],
                data['rev_primer'], data['product_size'],
                data['fwd_tm'], data['rev_tm'], data['amplicon_sequence'],
            ]]
            f.write(",".join(row) + "\n")
    print(f"\nValidated primer results saved → {output_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python design_amplicon.py <seed> <RANGES_CSV> "
              "<ONLY_GENOME_DIR> <BACKGROUND_DIR> <LENGTH> [MAX_LENGTH]")
        sys.exit(1)

    SEED            = int(sys.argv[1])
    RANGES_CSV      = sys.argv[2]
    ONLY_GENOME_DIR = sys.argv[3]
    BACKGROUND_DIR  = sys.argv[4]
    MIN_LENGTH      = int(sys.argv[5])
    MAX_LENGTH      = int(sys.argv[6]) if len(sys.argv) > 6 else MIN_LENGTH

    barcode_suffix = f"k{MIN_LENGTH}_{MAX_LENGTH}" if MAX_LENGTH != MIN_LENGTH else f"k{MIN_LENGTH}"
    ##barcode_suffix = "k60"
    CANDIDATE_DIR = f"./output/barcodes_{barcode_suffix}_amplicon_candidates/"
    OUTPUT_DIR    = f"./output/barcodes_{barcode_suffix}_amplicon_design/"

    min_hamming_dist = 3

    np.random.seed(SEED)
    random.seed(SEED)

    # ── Load candidate primers ────────────────────────────────────────────
    # Keys are (seq_id, region) tuples e.g. ("KL13_1", "Ori")
    all_candidates = load_all_candidate_primers(CANDIDATE_DIR)
    if not all_candidates:
        print(f"\nERROR: No candidate CSVs found in {CANDIDATE_DIR}.")
        sys.exit(1)
    print(f"Loaded candidates for {len(all_candidates)} (strain, region) pairs.")

    # ── Step 1: weighted combination selection ────────────────────────────
    candidate_counts, selection, failed_step1, failed_step1_causes = \
        select_optimal_combination_weighted(all_candidates, min_hamming_dist)

    # ── Step 2: cross-reactivity check ───────────────────────────────────
    print("Loading all background sequences for cross-reactivity check...")
    all_sequences = parse_all_records_from_dir(BACKGROUND_DIR)

    print("Checking cross-reactivity...", end="", flush=True)
    t_start = timer.perf_counter()
    problematic_pairs = check_cross_reactivity(
        selection, all_sequences,
        max_pcr_product=300, seed_len=4, max_mismatches=3,
    )
    print(f" done ({timer.perf_counter() - t_start:.2f}s)")

    failed_step2 = resolve_conflicts_greedy(problematic_pairs, candidate_counts)

    final_primers  = {k: data for k, data in selection.items()
                      if k not in failed_step2}
    current_count  = len(final_primers)
    print(f"{current_count} amplicons found on seed {SEED} (target: 36)")

    if current_count > 0:
        output_filename = f"validated_primers_{barcode_suffix}-seed{SEED}.csv"
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        output_path = os.path.join(OUTPUT_DIR, output_filename)
        with open(output_path, 'w', newline='') as f:
            header = ["Seq_ID", "Region", "Target_Barcode", "Forward_Primer",
                      "Reverse_Primer", "Product_Size", "Fwd_Tm", "Rev_Tm",
                      "Amplicon_Sequence"]
            f.write(",".join(header) + "\n")
            for (seq_id, region), data in sorted(final_primers.items()):
                row = [str(x) for x in [
                    seq_id, region, data['target_barcode'], data['fwd_primer'],
                    data['rev_primer'], data['product_size'],
                    data['fwd_tm'], data['rev_tm'], data['amplicon_sequence'],
                ]]
                f.write(",".join(row) + "\n")
        print(f"\n--- RUN COMPLETE ---")
        print(f"Seed {SEED}: {current_count} amplicons ({barcode_suffix})")
        print(f"Results → {output_path}")
    else:
        print(f"\n--- RUN COMPLETE ---")
        print(f"Seed {SEED}: no valid primer sets found.")

    if failed_step1 or failed_step2:
        print("\n--- FAILURE BREAKDOWN ---")
        if failed_step1:
            print(f"\nStep 1 (Candidate Selection): {len(failed_step1)} failed")
            for key, (h, d) in zip(failed_step1, failed_step1_causes):
                seq_id, region = key
                causes = []
                if h > 0: causes.append(f"barcode_hamming ({h})")
                if d > 0: causes.append(f"primer_dimer ({d})")
                if not causes: causes.append("no_candidates")
                print(f"  - {seq_id} {region}: {', '.join(causes)}")
        if failed_step2:
            print(f"\nStep 2 (Cross-Reactivity): {len(failed_step2)} failed")
            for seq_id, region in sorted(failed_step2):
                print(f"  - {seq_id} {region}")

    all_failed_ids = set(failed_step1) | failed_step2
