"""
rescue_amplicon.py — locked-set amplicon design
===============================================
One script for the two situations where part of the panel is already fixed:

  MODE "rescue"  Re-attempt targets that failed a previous design run, keeping
                 all previously validated primer pairs unchanged.

  MODE "extend"  Add newly acquired isolates/plasmids to an existing panel,
                 keeping every already-synthesized barcode/primer unchanged.

Both modes are the same operation: "design amplicons for target list T, given a
locked set L that must not change." The only differences are how T is chosen and
how hard the script audits L against sequences that were not in the original
reference set.

Pipeline
--------
  Step A  Audit the locked set against the expanded reference (extend mode).
          Three checks that the original pipeline never performs:
            A1. each locked barcode is still unique across old + new sequences
            A2. each locked primer pair still amplifies only its own target
            A3. locked_i x locked_j cross-combinations do not produce a product
                on the newly added templates
          Findings are reported, never silently fixed — the primers are already
          ordered, so the decision to drop a compromised amplicon is yours.
  Step 0  Generate fresh primer candidates for every target in T by sampling
          barcodes from the k-mer CSV (mandatory in extend mode: new isolates
          have no Stage-2 candidate files).
  Step 1  Pre-filter candidates against the locked set (Hamming + primer-dimer).
  Step 2  Joint weighted selection among the targets in T.
  Step 2.5 Re-verify each selected pair's own-target specificity against the
          full expanded reference (catches stale Stage-2 candidate files).
  Step 3  Cross-reactivity: new-vs-locked and new-vs-new.
  Step 4  Greedy conflict resolution that can only remove new targets.
          Optional --retry loop blacklists the failing candidate and re-selects.

Usage
-----
  # extend an existing panel with new isolates
  python rescue_amplicon.py --mode extend --comm comm3 --kmer 20 \
      --min-len 80 --max-len 120 --seed 1 \
      --locked ./output/comm3_k20_L80-120_amplicon_design/validated_primers_k20-seed1.csv \
      --new-ids ./new_isolates.txt

  # rescue targets that failed a previous run
  python rescue_amplicon.py --mode rescue --comm comm3 --kmer 20 \
      --min-len 80 --max-len 120 --seed 1 \
      --locked ./output/comm3_k20_L80-120_amplicon_design/validated_primers_k20-seed1.csv \
      --targets ./output/comm3_k20_L80-120_amplicon_design/failed_targets_k20-seed1.txt

  # audit only: is my existing panel still valid now that new plasmids exist?
  python rescue_amplicon.py --mode audit --comm comm3 --kmer 20 \
      --min-len 80 --max-len 120 --locked <validated.csv> --new-ids new_isolates.txt

Legacy positional form (equivalent to --mode rescue) is still accepted:
  python rescue_amplicon.py <seed> <COMM> <KMER> <MIN_LEN> <MAX_LEN> <locked.csv> <failed.txt>
"""

import argparse
import os
import random
import sys
import time as timer

import numpy as np
from Bio.Seq import Seq

from read_file_func import *
from design_amplicon import (
    calculate_hamming_distance,
    has_cross_dimer,
    select_optimal_combination_weighted,
    save_primer_results,
)
from individual_target_amplicon_candidates import design_candidate_primers
from sequence_alignment import is_primer_pair_specific


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_validated_primers(csv_path):
    """
    Reads back a validated_primers CSV (produced by save_primer_results) into
    the dict format used throughout the pipeline.
    """
    locked = {}
    with open(csv_path, 'r') as f:
        f.readline()  # header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            locked[parts[0]] = {
                'target_barcode':    parts[1],
                'fwd_primer':        parts[2],
                'rev_primer':        parts[3],
                'product_size':      int(parts[4]),
                'fwd_tm':            float(parts[5]),
                'rev_tm':            float(parts[6]),
                'amplicon_sequence': parts[7],
            }
    return locked


def load_id_list(txt_path):
    """Reads a plain-text file with one seq_id per line ('#' comments allowed)."""
    ids = []
    with open(txt_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line)
    return ids


# ---------------------------------------------------------------------------
# Barcode distance (length-safe)
# ---------------------------------------------------------------------------

_LENGTH_WARNED = set()

def barcodes_too_close(bc_a, bc_b, min_hamming_distance):
    """
    True if two barcodes are closer than min_hamming_distance.

    Hamming distance is only defined for equal-length strings. Barcodes of
    different lengths are treated as distinguishable (they are, by read length),
    but the mismatch is reported once so a k-mer-size mix-up cannot pass silently.
    """
    if len(bc_a) != len(bc_b):
        key = (len(bc_a), len(bc_b))
        if key not in _LENGTH_WARNED:
            _LENGTH_WARNED.add(key)
            print(f"  NOTE: comparing barcodes of length {len(bc_a)} and {len(bc_b)}. "
                  f"If this is unintended, the locked panel was built with a "
                  f"different --kmer than this run.")
        return False
    return calculate_hamming_distance(bc_a, bc_b) < min_hamming_distance


# ---------------------------------------------------------------------------
# Step A — audit the locked set against sequences it never saw
# ---------------------------------------------------------------------------

def count_motif_occurrences(sequences, motif):
    """Counts motif + reverse complement hits per sequence (linear, like kmer_generation)."""
    motif = motif.upper()
    rc = str(Seq(motif).reverse_complement()).upper()
    hits = {}
    for sid, seq in sequences.items():
        s = str(seq).upper()
        n = s.count(motif)
        if rc != motif:
            n += s.count(rc)
        if n:
            hits[sid] = n
    return hits


def audit_locked_set(locked_primers, all_sequences, new_sequences,
                     max_pcr_product=300, seed_len=12, max_mismatches=2,
                     check_locked_pairs=True):
    """
    Verifies that an already-synthesized panel is still valid after new templates
    were added to the reference set.

    Returns a dict:
        {'barcode_collisions': {seq_id: {other_seq_id: count}},
         'nonspecific':        [seq_id, ...],
         'cross_pairs':        [(seq_id_a, seq_id_b), ...]}
    """
    report = {'barcode_collisions': {}, 'nonspecific': [], 'cross_pairs': []}

    # A1 — barcode uniqueness across the expanded reference
    for seq_id, data in locked_primers.items():
        hits = count_motif_occurrences(all_sequences, data['target_barcode'])
        total = sum(hits.values())
        extra = {sid: n for sid, n in hits.items() if sid != seq_id}
        if total != 1 or extra:
            report['barcode_collisions'][seq_id] = hits

    # A2 — each locked pair still amplifies only its own target
    for seq_id, data in locked_primers.items():
        if not is_primer_pair_specific(
            data['fwd_primer'], data['rev_primer'], seq_id, all_sequences,
            max_pcr_product=max_pcr_product, seed_len=seed_len,
            max_mismatches=max_mismatches,
        ):
            report['nonspecific'].append(seq_id)

    # A3 — locked x locked cross-combinations on the NEW templates only
    # (locked-vs-locked was already cleared against the original reference, so
    #  restricting the search space to new templates keeps this cheap.)
    if check_locked_pairs and new_sequences:
        items = list(locked_primers.items())
        for i in range(len(items)):
            id_a, a = items[i]
            for j in range(i + 1, len(items)):
                id_b, b = items[j]
                combos = [
                    (a['fwd_primer'], b['fwd_primer']),
                    (a['fwd_primer'], b['rev_primer']),
                    (a['rev_primer'], b['fwd_primer']),
                    (a['rev_primer'], b['rev_primer']),
                ]
                for fwd, rev in combos:
                    if not is_primer_pair_specific(
                        fwd, rev, "", new_sequences,
                        max_pcr_product=max_pcr_product, seed_len=seed_len,
                        max_mismatches=max_mismatches,
                    ):
                        report['cross_pairs'].append((id_a, id_b))
                        break

    return report


def print_audit(report, n_locked, n_new):
    print(f"\n--- LOCKED-SET AUDIT ({n_locked} locked amplicons vs {n_new} new templates) ---")
    clean = True

    if report['barcode_collisions']:
        clean = False
        print(f"\n  [A1] {len(report['barcode_collisions'])} locked barcodes are no longer unique:")
        for sid, hits in sorted(report['barcode_collisions'].items()):
            where = ", ".join(f"{k} x{v}" for k, v in sorted(hits.items()))
            print(f"    - {sid}: found in {where}")

    if report['nonspecific']:
        clean = False
        print(f"\n  [A2] {len(report['nonspecific'])} locked primer pairs now amplify an off-target template:")
        for sid in sorted(report['nonspecific']):
            print(f"    - {sid}")

    if report['cross_pairs']:
        clean = False
        print(f"\n  [A3] {len(report['cross_pairs'])} locked primer pairs cross-amplify on the new templates:")
        for a, b in sorted(report['cross_pairs']):
            print(f"    - {a} x {b}")

    if clean:
        print("  All locked amplicons remain unique and specific. No action needed.")
    else:
        print("\n  These amplicons are already synthesized, so nothing was changed.")
        print("  Options: exclude the affected targets from analysis, re-design them")
        print("  (move their IDs into --targets), or drop the offending new isolate.")
    return clean


# ---------------------------------------------------------------------------
# Step 1 — pre-filter candidates against the locked set
# ---------------------------------------------------------------------------

def prefilter_candidates(locked_primers, new_candidates, min_hamming_distance):
    """
    Keeps only candidates compatible with ALL locked primers
    (barcode Hamming distance + primer-dimer).

    Returns (filtered_candidates, prefilter_failed_ids).
    """
    locked_list = list(locked_primers.values())
    filtered = {}
    prefilter_failed = []

    for seq_id, candidates in new_candidates.items():
        compatible = []
        for candidate in candidates:
            is_compatible = True
            for existing in locked_list:
                if barcodes_too_close(candidate['target_barcode'],
                                      existing['target_barcode'],
                                      min_hamming_distance):
                    is_compatible = False
                    break
                if has_cross_dimer(
                    candidate['fwd_primer'], candidate['rev_primer'],
                    existing['fwd_primer'], existing['rev_primer'],
                    max_dimer_dg=-9000.0, max_3prime_dg=-9000.0,
                ):
                    is_compatible = False
                    break
            if is_compatible:
                compatible.append(candidate)

        if compatible:
            filtered[seq_id] = compatible
        else:
            prefilter_failed.append(seq_id)

    return filtered, prefilter_failed


# ---------------------------------------------------------------------------
# Step 3 — cross-reactivity (skips locked-vs-locked; that is Step A3's job)
# ---------------------------------------------------------------------------

def check_cross_reactivity_new(locked_primers, new_primers, all_sequences,
                               max_pcr_product=300, seed_len=12, max_mismatches=2):
    """Checks new-vs-locked and new-vs-new primer combinations."""
    problematic_pairs = set()
    locked_list = list(locked_primers.items())
    new_list = list(new_primers.items())

    def _check_pair(primer_A, primer_B):
        combos = [
            (primer_A['fwd_primer'], primer_B['fwd_primer']),
            (primer_A['fwd_primer'], primer_B['rev_primer']),
            (primer_A['rev_primer'], primer_B['fwd_primer']),
            (primer_A['rev_primer'], primer_B['rev_primer']),
        ]
        for fwd, rev in combos:
            if not is_primer_pair_specific(
                fwd, rev, "", all_sequences,
                max_pcr_product=max_pcr_product, seed_len=seed_len,
                max_mismatches=max_mismatches,
            ):
                return True
        return False

    for n_id, n_data in new_list:
        for l_id, l_data in locked_list:
            if _check_pair(n_data, l_data):
                problematic_pairs.add((n_id, l_id))

    for i in range(len(new_list)):
        for j in range(i + 1, len(new_list)):
            a_id, a_data = new_list[i]
            b_id, b_data = new_list[j]
            if _check_pair(a_data, b_data):
                problematic_pairs.add((a_id, b_id))

    return problematic_pairs


def resolve_conflicts_locked_safe(problematic_pairs, locked_ids, candidate_counts=None):
    """
    Greedy conflict resolution that never removes a locked target.
    Tie-break: prefer keeping the target with fewer candidates available.
    """
    if not problematic_pairs:
        return set()

    graph = {}
    for a, b in problematic_pairs:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)

    alive = set(graph.keys())
    failed_ids = set()
    counts = candidate_counts or {}

    def degree(n):
        return len([nbr for nbr in graph.get(n, set()) if nbr in alive])

    while True:
        candidates = [n for n in alive if degree(n) > 0 and n not in locked_ids]
        if not candidates:
            break
        pick = max(candidates, key=lambda n: (degree(n), -counts.get(n, 0)))
        alive.remove(pick)
        failed_ids.add(pick)

    remaining = [n for n in alive if degree(n) > 0]
    if remaining:
        print(f"  WARNING: {len(remaining)} locked targets still have unresolved "
              f"cross-reactivity conflicts (they were kept anyway): {sorted(remaining)}")

    return failed_ids


# ---------------------------------------------------------------------------
# Steps 2-4 — selection + cross-reactivity, single pass or retry loop
# ---------------------------------------------------------------------------

def select_and_validate(filtered_candidates, locked_primers, sequences,
                        min_hamming_dist, max_rounds=1,
                        max_pcr_product=300, seed_len=12, max_mismatches=2,
                        verify_own_specificity=True):
    """
    Runs joint selection followed by the cross-reactivity check.

    max_rounds=1  -> single pass (failed targets are dropped)
    max_rounds>1  -> retry loop (the failing candidate is blacklisted and the
                     target competes again in the next round)

    Returns (accepted, selection_failed, exhausted, still_pending).
    """
    accepted = {}
    current = {sid: list(c) for sid, c in filtered_candidates.items()}
    selection_failed = set()
    exhausted = set()
    round_num = 0

    while current and round_num < max_rounds:
        round_num += 1
        current_locked = {**locked_primers, **accepted}
        if max_rounds > 1:
            print(f"\n  === Round {round_num}: {len(current)} targets, "
                  f"{sum(len(v) for v in current.values())} candidates ===")

        t0 = timer.perf_counter()
        candidate_counts, selection, failed_sel, _ = \
            select_optimal_combination_weighted(current, min_hamming_dist)
        print(f"  Selection: {len(selection)} selected, {len(failed_sel)} failed "
              f"({timer.perf_counter() - t0:.2f}s)")
        selection_failed.update(failed_sel)

        if not selection:
            break

        # Step 2.5 — own-target specificity against the full expanded reference.
        if verify_own_specificity:
            stale = []
            for sid, data in list(selection.items()):
                if not is_primer_pair_specific(
                    data['fwd_primer'], data['rev_primer'], sid, sequences,
                    max_pcr_product=max_pcr_product, seed_len=seed_len,
                    max_mismatches=max_mismatches,
                ):
                    stale.append(sid)
            if stale:
                print(f"  {len(stale)} selected pairs amplify an off-target template "
                      f"and were dropped: {sorted(stale)}")
                for sid in stale:
                    key = (selection[sid]['fwd_primer'], selection[sid]['rev_primer'])
                    if sid in current:
                        current[sid] = [c for c in current[sid]
                                        if (c['fwd_primer'], c['rev_primer']) != key]
                        if not current[sid]:
                            del current[sid]
                            exhausted.add(sid)
                    del selection[sid]

        if not selection:
            continue

        print("  Cross-reactivity...", end="", flush=True)
        t0 = timer.perf_counter()
        problematic_pairs = check_cross_reactivity_new(
            current_locked, selection, sequences,
            max_pcr_product=max_pcr_product, seed_len=seed_len,
            max_mismatches=max_mismatches,
        )
        print(f" done ({timer.perf_counter() - t0:.2f}s)")

        if not problematic_pairs:
            accepted.update(selection)
            for sid in selection:
                current.pop(sid, None)
                selection_failed.discard(sid)
            print(f"  All {len(selection)} passed. Total accepted: {len(accepted)}")
            break

        failed_xr = resolve_conflicts_locked_safe(
            problematic_pairs, set(current_locked.keys()), candidate_counts
        )
        survived = {pid: d for pid, d in selection.items() if pid not in failed_xr}
        accepted.update(survived)
        for sid in survived:
            current.pop(sid, None)
            selection_failed.discard(sid)
        print(f"  {len(failed_xr)} failed cross-reactivity, {len(survived)} survived.")

        # Blacklist the specific candidate that failed, then retry if allowed.
        retryable = 0
        for sid in failed_xr:
            key = (selection[sid]['fwd_primer'], selection[sid]['rev_primer'])
            if sid in current:
                current[sid] = [c for c in current[sid]
                                if (c['fwd_primer'], c['rev_primer']) != key]
                if current[sid]:
                    retryable += 1
                    selection_failed.discard(sid)
                else:
                    del current[sid]
                    exhausted.add(sid)

        if max_rounds == 1:
            exhausted.update(failed_xr - set(accepted.keys()))
            break
        if retryable == 0 and not failed_sel:
            print("  No retryable targets. Stopping.")
            break

    still_pending = set(current.keys())
    return accepted, (selection_failed - set(accepted.keys())), exhausted, still_pending


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv):
    # Legacy positional form: seed COMM KMER MIN MAX locked.csv failed.txt
    if len(argv) == 8 and not argv[1].startswith("-"):
        argv = [argv[0],
                "--mode", "rescue",
                "--seed", argv[1], "--comm", argv[2], "--kmer", argv[3],
                "--min-len", argv[4], "--max-len", argv[5],
                "--locked", argv[6], "--targets", argv[7]]

    p = argparse.ArgumentParser(description="Locked-set amplicon design (rescue / extend / audit).")
    p.add_argument("--mode", choices=["rescue", "extend", "audit"], default="rescue")
    p.add_argument("--comm", required=True, help="community name; targets live in ./targets/<comm>")
    p.add_argument("--kmer", type=int, required=True, help="barcode length (must match the locked panel)")
    p.add_argument("--min-len", type=int, required=True)
    p.add_argument("--max-len", type=int, required=True)
    p.add_argument("--locked", required=True, help="validated_primers CSV to keep fixed")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--targets", default=None,
                   help="file of seq_ids to design (default in extend mode: every "
                        "sequence in ./targets/<comm> that is not in --locked)")
    p.add_argument("--new-ids", default=None,
                   help="file of newly added seq_ids; scopes the locked-set audit")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--min-hamming", type=int, default=2)
    p.add_argument("--num-barcodes", type=int, default=100,
                   help="barcodes sampled per target during candidate generation")
    p.add_argument("--retry", action="store_true", help="enable the multi-round retry loop")
    p.add_argument("--rounds", type=int, default=10, help="max rounds when --retry is set")
    p.add_argument("--use-existing-candidates", action="store_true",
                   help="also load Stage-2 candidate CSVs (only safe if they were "
                        "generated against the current, complete reference set)")
    p.add_argument("--no-fresh-candidates", action="store_true",
                   help="do not run Primer3; use Stage-2 candidate CSVs only")
    p.add_argument("--skip-locked-audit", action="store_true")
    p.add_argument("--skip-locked-pair-audit", action="store_true",
                   help="skip check A3 (locked x locked on new templates); it is the slow one")
    p.add_argument("--max-pcr-product", type=int, default=300)
    p.add_argument("--seed-len", type=int, default=12)
    p.add_argument("--max-mismatches", type=int, default=2)
    return p.parse_args(argv[1:])


def main(argv):
    args = parse_args(argv)

    SEQ_DIR = f"./targets/{args.comm}"
    BARCODE_CSV = f"./output/KMERS/{args.comm}_k{args.kmer}.csv"
    CANDIDATE_DIR = (f"./output/KMERS/{args.comm}_k{args.kmer}"
                     f"_L{args.min_len}-{args.max_len}_amplicon_candidates/")
    OUTPUT_DIR = args.out_dir or (f"./output/{args.comm}_k{args.kmer}"
                                  f"_L{args.min_len}-{args.max_len}_amplicon_{args.mode}/")

    np.random.seed(args.seed)
    random.seed(args.seed)

    primer_config = {
        'flank_size': 300,
        'primer3_globals': {
            'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0, 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0, 'PRIMER_MAX_GC': 60.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[args.min_len, args.max_len]],
        }
    }

    # --- Load reference ---
    print(f"Loading sequences from {SEQ_DIR} ...")
    sequences = parse_sequences(find_fasta_files(SEQ_DIR))
    if not sequences:
        print("\nERROR: No sequences were parsed. Check SEQ_DIR.")
        return 1
    print(f"  {len(sequences)} sequences in the reference set.")

    print(f"Loading locked primers from {args.locked} ...")
    locked_primers = load_validated_primers(args.locked)
    locked_ids = set(locked_primers.keys())
    print(f"  {len(locked_primers)} locked amplicons.")

    missing = locked_ids - set(sequences.keys())
    if missing:
        print(f"  WARNING: {len(missing)} locked seq_ids have no FASTA in {SEQ_DIR}: "
              f"{sorted(missing)}")

    bad_len = {sid for sid, d in locked_primers.items()
               if len(d['target_barcode']) != args.kmer}
    if bad_len:
        print(f"  WARNING: {len(bad_len)} locked barcodes are not {args.kmer} bp. "
              f"Run with --kmer set to the value used for the locked panel.")

    # --- Which sequences are new? ---
    if args.new_ids:
        new_ids = [i for i in load_id_list(args.new_ids) if i in sequences]
    else:
        new_ids = [i for i in sequences if i not in locked_ids]
    new_sequences = {i: sequences[i] for i in new_ids}

    # --- Step A: audit the locked set ---
    if not args.skip_locked_audit:
        print(f"\nStep A: Auditing locked set against the expanded reference ...")
        t0 = timer.perf_counter()
        report = audit_locked_set(
            locked_primers, sequences, new_sequences,
            max_pcr_product=args.max_pcr_product, seed_len=args.seed_len,
            max_mismatches=args.max_mismatches,
            check_locked_pairs=not args.skip_locked_pair_audit,
        )
        print_audit(report, len(locked_primers), len(new_sequences))
        print(f"  (audit took {timer.perf_counter() - t0:.1f}s)")

        os.makedirs(OUTPUT_DIR, exist_ok=True)
        audit_path = os.path.join(OUTPUT_DIR, f"locked_audit_k{args.kmer}.txt")
        with open(audit_path, 'w') as f:
            for sid, hits in sorted(report['barcode_collisions'].items()):
                f.write(f"barcode_not_unique\t{sid}\t"
                        f"{';'.join(f'{k}:{v}' for k, v in sorted(hits.items()))}\n")
            for sid in sorted(report['nonspecific']):
                f.write(f"pair_not_specific\t{sid}\n")
            for a, b in sorted(report['cross_pairs']):
                f.write(f"locked_cross_pair\t{a}\t{b}\n")
        print(f"  Audit written to: {audit_path}")

    if args.mode == "audit":
        return 0

    # --- Target list ---
    if args.targets:
        target_ids = load_id_list(args.targets)
    elif args.mode == "extend":
        target_ids = sorted(set(sequences.keys()) - locked_ids)
    else:
        print("\nERROR: --targets is required in rescue mode.")
        return 1

    target_ids = [t for t in target_ids if t not in locked_ids]
    print(f"\n{len(target_ids)} targets to design "
          f"({'extend' if args.mode == 'extend' else 'rescue'} mode).")
    if not target_ids:
        print("Nothing to do.")
        return 0

    # --- Step 0: candidate generation ---
    print(f"\nStep 0: Building candidate pools ...")
    t0 = timer.perf_counter()
    existing = load_all_candidate_primers(CANDIDATE_DIR) if args.use_existing_candidates else {}
    new_candidates = {}
    no_candidates = []
    exhausted_pool = []   # targets whose entire barcode space is sampled every seed

    for seq_id in target_ids:
        old = existing.get(seq_id, [])
        fresh = []
        n_barcodes = 0
        if not args.no_fresh_candidates:
            barcodes = load_unique_barcodes(BARCODE_CSV, seq_id)
            n_barcodes = len(barcodes)
            if barcodes:
                if n_barcodes <= args.num_barcodes:
                    exhausted_pool.append(seq_id)
                fresh = design_candidate_primers(
                    sequences, seq_id, barcodes, primer_config,
                    num_of_barcodes_checked=args.num_barcodes,
                )
            else:
                print(f"  {seq_id}: no unique barcodes in {BARCODE_CSV}")
        combined = old + fresh
        if combined:
            new_candidates[seq_id] = combined
            flag = " [pool exhausted]" if seq_id in exhausted_pool else ""
            print(f"  {seq_id}: {len(old)} existing + {len(fresh)} fresh = {len(combined)} "
                  f"(from {n_barcodes} unique barcodes){flag}")
        else:
            no_candidates.append(seq_id)
            print(f"  {seq_id}: no candidates (from {n_barcodes} unique barcodes)")
    print(f"  Candidate generation done ({timer.perf_counter() - t0:.1f}s)")

    if exhausted_pool:
        print(f"\n  NOTE: {len(exhausted_pool)} targets have <= --num-barcodes "
              f"({args.num_barcodes}) unique barcodes, so every seed samples their entire "
              f"barcode space and produces the identical candidate pool:")
        print(f"    {sorted(exhausted_pool)}")
        print("  Re-running with more seeds cannot help these. Raise --kmer to create more "
              "unique barcodes, widen the amplicon size range, or relax the Primer3 limits.")

    if not new_candidates:
        print("\nNo candidates available. Nothing to add.")
        return 0

    # --- Step 1: pre-filter against locked ---
    print(f"\nStep 1: Pre-filtering against {len(locked_primers)} locked primers ...")
    t0 = timer.perf_counter()
    filtered_candidates, prefilter_failed = prefilter_candidates(
        locked_primers, new_candidates, args.min_hamming
    )
    before = sum(len(v) for v in new_candidates.values())
    after = sum(len(v) for v in filtered_candidates.values())
    print(f"  {before} -> {after} candidates compatible with the locked set "
          f"({timer.perf_counter() - t0:.1f}s)")
    print(f"  {len(filtered_candidates)} targets retained, {len(prefilter_failed)} eliminated")

    # --- Steps 2-4 ---
    if filtered_candidates:
        print(f"\nSteps 2-4: Selection and cross-reactivity "
              f"({'retry loop' if args.retry else 'single pass'}) ...")
        accepted, selection_failed, exhausted, still_pending = select_and_validate(
            filtered_candidates, locked_primers, sequences,
            args.min_hamming,
            max_rounds=args.rounds if args.retry else 1,
            max_pcr_product=args.max_pcr_product, seed_len=args.seed_len,
            max_mismatches=args.max_mismatches,
        )
    else:
        accepted, selection_failed, exhausted, still_pending = {}, set(), set(), set()

    # --- Merge and save ---
    merged = {**locked_primers, **accepted}
    still_failed = (set(no_candidates) | set(prefilter_failed) |
                    selection_failed | exhausted | still_pending) - set(accepted.keys())

    print(f"\n--- SUMMARY ({args.mode}, seed {args.seed}) ---")
    print(f"  Locked (unchanged):  {len(locked_primers)}")
    print(f"  Newly designed:      {len(accepted)}")
    print(f"  Final panel size:    {len(merged)}")
    print(f"  Unresolved targets:  {len(still_failed)}")

    prefix = "extended" if args.mode == "extend" else "rescued"
    output_filename = f"{prefix}_primers_k{args.kmer}-seed{args.seed}.csv"
    save_primer_results(OUTPUT_DIR, merged, filename=output_filename)

    if accepted:
        save_primer_results(OUTPUT_DIR, accepted,
                            filename=f"{prefix}_primers_NEW_ONLY_k{args.kmer}-seed{args.seed}.csv")

    if still_failed:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        failed_path = os.path.join(OUTPUT_DIR, f"still_failed_k{args.kmer}-seed{args.seed}.txt")
        with open(failed_path, 'w') as f:
            for sid in sorted(still_failed):
                f.write(sid + "\n")
        print(f"  Unresolved target IDs: {failed_path}")

        print("\n--- FAILURE BREAKDOWN ---")
        if no_candidates:
            print(f"  No candidates at all ({len(no_candidates)}): {sorted(no_candidates)}")
        if prefilter_failed:
            print(f"  Incompatible with locked set ({len(prefilter_failed)}): {sorted(prefilter_failed)}")
        sel = sorted(selection_failed - exhausted)
        if sel:
            print(f"  Mutually incompatible among new targets ({len(sel)}): {sel}")
        if exhausted:
            print(f"  Cross-reactivity, candidates exhausted ({len(exhausted)}): {sorted(exhausted)}")
        if still_pending - exhausted - selection_failed:
            print(f"  Round limit reached: {sorted(still_pending - exhausted - selection_failed)}")

    print("\n--- COMPLETE ---")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
