"""
Add amplicons to an existing panel, keeping already-designed primers fixed.

Usage
-----
  # add every isolate that is not already in the locked panel
  python add_amplicons.py \
      --locked ./output/barcodes_k20_amplicon_design/validated_primers_k20-seed7.csv \
      --only-genome ./only_genome \
      --background ./genome_and_plasmids_within_host \
      --length 20 --seed 1

  # audit only: did the new isolates break the existing panel?
  python add_amplicons.py --mode audit --locked <csv> \
      --only-genome ./only_genome --background ./genome_and_plasmids_within_host \
      --length 20

  # design only specific targets
  printf 'KL30_1,Ori\nKL30_1,Ter\nKL13_1,Ter\n' > targets.txt
  python add_amplicons.py --locked <csv> --targets targets.txt ...
"""

import argparse
import csv
import os
import random
import re
import sys
import time as timer

# Ensure the script's own directory is on the path so local modules are found
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

from read_file_func import (
    find_fasta_files,
    parse_sequences,
    load_unique_barcodes,
    load_all_candidate_primers,
)
from design_amplicon import (
    calculate_hamming_distance,
    has_cross_dimer,
    select_optimal_combination_weighted,
)
from individual_target_amplicon_candidates import design_candidate_primers
from sequence_alignment import is_primer_pair_specific


REGIONS = ("Ori", "Ter")


def load_locked_panel(csv_path):
    """
    Read a validated_primers CSV into {(seq_id, region): candidate_dict},
    the same key shape Stage 3 uses. Tolerates a legacy file with no Region
    column by assigning region ''.
    """
    locked = {}
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        has_region = "Region" in (reader.fieldnames or [])
        for row in reader:
            seq_id = row["Seq_ID"].strip()
            region = row["Region"].strip() if has_region else ""
            locked[(seq_id, region)] = {
                "target_barcode":    row["Target_Barcode"],
                "fwd_primer":        row["Forward_Primer"],
                "rev_primer":        row["Reverse_Primer"],
                "product_size":      int(row["Product_Size"]),
                "fwd_tm":            float(row["Fwd_Tm"]),
                "rev_tm":            float(row["Rev_Tm"]),
                "amplicon_sequence": row["Amplicon_Sequence"],
            }
    if not has_region:
        print("  WARNING: locked CSV has no Region column; all entries assigned "
              "region ''. Regenerate it with the current design_amplicon.py.")
    return locked


def save_panel(output_dir, primers, filename):
    """Write {(seq_id, region): data} in the Stage 3 column layout."""
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    with open(output_path, "w", newline="") as f:
        header = ["Seq_ID", "Region", "Target_Barcode", "Forward_Primer",
                  "Reverse_Primer", "Product_Size", "Fwd_Tm", "Rev_Tm",
                  "Amplicon_Sequence"]
        f.write(",".join(header) + "\n")
        for (seq_id, region), data in sorted(primers.items()):
            row = [str(x) for x in [
                seq_id, region, data["target_barcode"], data["fwd_primer"],
                data["rev_primer"], data["product_size"],
                data["fwd_tm"], data["rev_tm"], data["amplicon_sequence"],
            ]]
            f.write(",".join(row) + "\n")
    print(f"  Saved {len(primers)} amplicons -> {output_path}")
    return output_path


def load_targets_file(txt_path):
    """
    One target per line. Either 'SEQ_ID,Region' or bare 'SEQ_ID'
    (which expands to both regions).
    """
    targets = []
    with open(txt_path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if "," in line:
                seq_id, region = (p.strip() for p in line.split(",", 1))
                targets.append((seq_id, region))
            else:
                targets.extend((line, r) for r in REGIONS)
    return targets


def parse_background_by_file(directory):
    all_records = {}
    stem_to_ids = {}
    dir_path = Path(directory)
    fasta_files = sorted(
        [f for f in dir_path.glob("*.fasta") if not f.name.startswith("._")] +
        [f for f in dir_path.glob("*.fa") if not f.name.startswith("._")]
    )
    for fasta_file in fasta_files:
        ids = []
        for rec in SeqIO.parse(str(fasta_file), "fasta"):
            all_records[rec.id] = str(rec.seq).upper()
            ids.append(rec.id)
        stem_to_ids[fasta_file.stem] = ids
    print(f"Loaded {len(all_records)} records from {len(fasta_files)} files "
          f"in '{directory}'.")
    return all_records, stem_to_ids


def isolate_stem_for(seq_id, stem_to_ids):
    """
    Find the genome_and_plasmids file that contains a given chromosome record.
    """
    for stem, ids in stem_to_ids.items():
        if seq_id in ids:
            return stem
    return re.sub(r"_\d+$", "", seq_id)


### Barcode Distance
def barcodes_too_close(bc_a, bc_b, min_hamming_distance):
    """
    True if two barcodes are closer than min_hamming_distance.
    """
    dist = calculate_hamming_distance(bc_a, bc_b)
    return dist is not None and dist < min_hamming_distance


### Locked panel audit
def count_motif_occurrences(sequences, motif):
    """Count motif + reverse complement hits per record."""
    motif = motif.upper()
    rc = str(Seq(motif).reverse_complement()).upper()
    hits = {}
    for rec_id, seq in sequences.items():
        s = str(seq).upper()
        n = s.count(motif)
        if rc != motif:
            n += s.count(rc)
        if n:
            hits[rec_id] = n
    return hits


def audit_locked_panel(locked, all_records, new_record_ids,
                       max_pcr_product=300, seed_len=9, max_mismatches=2,
                       check_locked_pairs=True):
    """
    Returns {'unmapped': [...], 'barcode_collisions': {...},
             'nonspecific': [...], 'cross_pairs': [...]}
    """
    report = {"unmapped": [], "barcode_collisions": {}, "nonspecific": [],
              "cross_pairs": []}

    # A0 - every locked Seq_ID must resolve to a record
    for (seq_id, region) in locked:
        if seq_id not in all_records:
            report["unmapped"].append((seq_id, region))

    # A1 — barcode still occurs exactly once
    for (seq_id, region), data in locked.items():
        hits = count_motif_occurrences(all_records, data["target_barcode"])
        total = sum(hits.values())
        elsewhere = {r: n for r, n in hits.items() if r != seq_id}
        if total != 1 or elsewhere:
            report["barcode_collisions"][(seq_id, region)] = hits

    # A2 — pair still amplifies only its own target
    for (seq_id, region), data in locked.items():
        if seq_id not in all_records:
            continue 
        if not is_primer_pair_specific(
            data["fwd_primer"], data["rev_primer"], seq_id, all_records,
            max_pcr_product=max_pcr_product, seed_len=seed_len,
            max_mismatches=max_mismatches,
        ):
            report["nonspecific"].append((seq_id, region))

    # A3 — locked x locked, searched on the new records only.
    if check_locked_pairs and new_record_ids:
        new_records = {r: all_records[r] for r in new_record_ids if r in all_records}
        items = list(locked.items())
        for i in range(len(items)):
            key_a, a = items[i]
            for j in range(i + 1, len(items)):
                key_b, b = items[j]
                combos = [
                    (a["fwd_primer"], b["fwd_primer"]),
                    (a["fwd_primer"], b["rev_primer"]),
                    (a["rev_primer"], b["fwd_primer"]),
                    (a["rev_primer"], b["rev_primer"]),
                ]
                for fwd, rev in combos:
                    if not is_primer_pair_specific(
                        fwd, rev, "", new_records,
                        max_pcr_product=max_pcr_product, seed_len=seed_len,
                        max_mismatches=max_mismatches,
                    ):
                        report["cross_pairs"].append((key_a, key_b))
                        break

    return report


def _fmt(key):
    seq_id, region = key
    return f"{seq_id} {region}".strip()


def print_audit(report, n_locked, n_new):
    print(f"\n--- LOCKED-PANEL AUDIT ({n_locked} locked amplicons, "
          f"{n_new} newly added records) ---")
    clean = True

    if report["unmapped"]:
        clean = False
        print(f"\n  [A0] {len(report['unmapped'])} locked Seq_IDs have no matching "
              f"record in the background directory:")
        for k in sorted(report["unmapped"]):
            print(f"    - {_fmt(k)}")
        print("    Checks A1/A2 cannot exclude their own amplicon, so their "
              "results below are unreliable. Confirm that Seq_ID matches the "
              "chromosome record id (e.g. 'KL13_1').")

    if report["barcode_collisions"]:
        clean = False
        print(f"\n  [A1] {len(report['barcode_collisions'])} locked barcodes are "
              f"no longer unique:")
        for k, hits in sorted(report["barcode_collisions"].items()):
            where = ", ".join(f"{r} x{n}" for r, n in sorted(hits.items()))
            print(f"    - {_fmt(k)}: found in {where}")

    if report["nonspecific"]:
        clean = False
        print(f"\n  [A2] {len(report['nonspecific'])} locked primer pairs now "
              f"amplify an off-target record:")
        for k in sorted(report["nonspecific"]):
            print(f"    - {_fmt(k)}")

    if report["cross_pairs"]:
        clean = False
        print(f"\n  [A3] {len(report['cross_pairs'])} locked primer pairs "
              f"cross-amplify on a new record:")
        for a, b in sorted(report["cross_pairs"]):
            print(f"    - {_fmt(a)}  x  {_fmt(b)}")

    if clean:
        print("  All locked amplicons remain unique and specific. No action needed.")
    else:
        print("\n  Nothing was changed — these primers already exist. Options:")
        print("    - exclude the affected amplicon from quantification")
        print("    - add its target to --targets so a replacement is designed")
        print("    - leave the offending new isolate out of the community")
    return clean


def write_audit(report, output_dir, suffix):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"locked_audit_{suffix}.txt")
    with open(path, "w") as f:
        for k in sorted(report["unmapped"]):
            f.write(f"seq_id_not_in_background\t{k[0]}\t{k[1]}\n")
        for k, hits in sorted(report["barcode_collisions"].items()):
            f.write(f"barcode_not_unique\t{k[0]}\t{k[1]}\t"
                    f"{';'.join(f'{r}:{n}' for r, n in sorted(hits.items()))}\n")
        for k in sorted(report["nonspecific"]):
            f.write(f"pair_not_specific\t{k[0]}\t{k[1]}\n")
        for a, b in sorted(report["cross_pairs"]):
            f.write(f"locked_cross_pair\t{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\n")
    print(f"  Audit written to: {path}")
    return path



### Pre-filter against the locked panel

def prefilter_candidates(locked, new_candidates, min_hamming_distance):
    """Keep only candidates compatible with every locked amplicon."""
    locked_list = list(locked.values())
    filtered = {}
    prefilter_failed = []

    for key, candidates in new_candidates.items():
        compatible = []
        for cand in candidates:
            ok = True
            for existing in locked_list:
                if barcodes_too_close(cand["target_barcode"],
                                      existing["target_barcode"],
                                      min_hamming_distance):
                    ok = False
                    break
                if has_cross_dimer(cand["fwd_primer"], cand["rev_primer"],
                                   existing["fwd_primer"], existing["rev_primer"]):
                    ok = False
                    break
            if ok:
                compatible.append(cand)
        if compatible:
            filtered[key] = compatible
        else:
            prefilter_failed.append(key)

    return filtered, prefilter_failed


### Cross-reactivity check

def check_cross_reactivity_new(locked, new_primers, all_records,
                               max_pcr_product=300, seed_len=9, max_mismatches=2):
    problematic_pairs = set()
    locked_list = list(locked.items())
    new_list = list(new_primers.items())

    def _conflicts(a, b):
        combos = [
            (a["fwd_primer"], b["fwd_primer"]),
            (a["fwd_primer"], b["rev_primer"]),
            (a["rev_primer"], b["fwd_primer"]),
            (a["rev_primer"], b["rev_primer"]),
        ]
        for fwd, rev in combos:
            if not is_primer_pair_specific(
                fwd, rev, "", all_records,
                max_pcr_product=max_pcr_product, seed_len=seed_len,
                max_mismatches=max_mismatches,
            ):
                return True
        return False

    for n_key, n_data in new_list:
        for l_key, l_data in locked_list:
            if _conflicts(n_data, l_data):
                problematic_pairs.add((n_key, l_key))

    for i in range(len(new_list)):
        for j in range(i + 1, len(new_list)):
            a_key, a_data = new_list[i]
            b_key, b_data = new_list[j]
            if _conflicts(a_data, b_data):
                problematic_pairs.add((a_key, b_key))

    return problematic_pairs


def resolve_conflicts_locked_safe(problematic_pairs, locked_keys, candidate_counts=None):
    """Greedy resolution that can only remove new targets."""
    if not problematic_pairs:
        return set()

    graph = {}
    for a, b in problematic_pairs:
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)

    alive = set(graph.keys())
    failed = set()
    counts = candidate_counts or {}

    def degree(n):
        return sum(1 for nbr in graph.get(n, set()) if nbr in alive)

    while True:
        removable = [n for n in alive if degree(n) > 0 and n not in locked_keys]
        if not removable:
            break
        pick = max(removable, key=lambda n: (degree(n), -counts.get(n, 0)))
        alive.remove(pick)
        failed.add(pick)

    stuck = [n for n in alive if degree(n) > 0]
    if stuck:
        print(f"  WARNING: {len(stuck)} locked amplicons conflict with each other "
              f"and were kept anyway: {[_fmt(k) for k in sorted(stuck)]}")
    return failed


def select_and_validate(filtered_candidates, locked, all_records, min_hamming_dist,
                        max_rounds=1, max_pcr_product=300, seed_len=9,
                        max_mismatches=2):
    """
    Joint selection followed by cross-reactivity.

    Returns (accepted, selection_failed, exhausted, still_pending).
    """
    accepted = {}
    current = {k: list(v) for k, v in filtered_candidates.items()}
    selection_failed = set()
    exhausted = set()
    round_num = 0

    while current and round_num < max_rounds:
        round_num += 1
        current_locked = {**locked, **accepted}
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

        # Own-target specificity against the current reference
        stale = []
        for key, data in list(selection.items()):
            seq_id, _ = key
            if not is_primer_pair_specific(
                data["fwd_primer"], data["rev_primer"], seq_id, all_records,
                max_pcr_product=max_pcr_product, seed_len=seed_len,
                max_mismatches=max_mismatches,
            ):
                stale.append(key)
        if stale:
            print(f"  {len(stale)} selected pairs amplify an off-target record "
                  f"and were dropped: {[_fmt(k) for k in stale]}")
            for key in stale:
                bl = (selection[key]["fwd_primer"], selection[key]["rev_primer"])
                if key in current:
                    current[key] = [c for c in current[key]
                                    if (c["fwd_primer"], c["rev_primer"]) != bl]
                    if not current[key]:
                        del current[key]
                        exhausted.add(key)
                del selection[key]

        if not selection:
            continue

        print("  Cross-reactivity...", end="", flush=True)
        t0 = timer.perf_counter()
        problematic_pairs = check_cross_reactivity_new(
            current_locked, selection, all_records,
            max_pcr_product=max_pcr_product, seed_len=seed_len,
            max_mismatches=max_mismatches,
        )
        print(f" done ({timer.perf_counter() - t0:.2f}s)")

        if not problematic_pairs:
            accepted.update(selection)
            for key in selection:
                current.pop(key, None)
                selection_failed.discard(key)
            print(f"  All {len(selection)} passed. Total accepted: {len(accepted)}")
            break

        failed_xr = resolve_conflicts_locked_safe(
            problematic_pairs, set(current_locked.keys()), candidate_counts
        )
        survived = {k: v for k, v in selection.items() if k not in failed_xr}
        accepted.update(survived)
        for key in survived:
            current.pop(key, None)
            selection_failed.discard(key)
        print(f"  {len(failed_xr)} failed cross-reactivity, {len(survived)} survived.")

        retryable = 0
        for key in failed_xr:
            bl = (selection[key]["fwd_primer"], selection[key]["rev_primer"])
            if key in current:
                current[key] = [c for c in current[key]
                                if (c["fwd_primer"], c["rev_primer"]) != bl]
                if current[key]:
                    retryable += 1
                    selection_failed.discard(key)
                else:
                    del current[key]
                    exhausted.add(key)

        if max_rounds == 1:
            exhausted.update(failed_xr - set(accepted.keys()))
            break
        if retryable == 0 and not failed_sel:
            print("  No retryable targets. Stopping.")
            break

    return accepted, (selection_failed - set(accepted.keys())), exhausted, set(current.keys())


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Add amplicons to an existing panel without changing it.")
    p.add_argument("--locked", required=True,
                   help="validated_primers CSV whose amplicons must not change")
    p.add_argument("--only-genome", default="./only_genome",
                   help="chromosome-only FASTA directory (primer design templates)")
    p.add_argument("--background", default="./genome_and_plasmids_within_host",
                   help="genome+plasmid FASTA directory (specificity reference)")
    p.add_argument("--length", type=int, required=True,
                   help="barcode MIN_LENGTH used in Stage 1")
    p.add_argument("--max-length", type=int, default=None,
                   help="barcode MAX_LENGTH used in Stage 1 (defaults to --length)")
    p.add_argument("--mode", choices=["add", "audit"], default="add")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--targets", default=None,
                   help="file of targets to design: 'SEQ_ID,Region' or bare SEQ_ID "
                        "per line. Default: every only_genome isolate x {Ori,Ter} "
                        "that is not already locked")
    p.add_argument("--unlock", default=None,
                   help="file of locked targets to release and re-design "
                        "(same format as --targets). Use this to replace an "
                        "amplicon the audit flagged as broken.")
    p.add_argument("--new-isolates", default=None,
                   help="file of newly added Seq_IDs, one per line; scopes audit A3. "
                        "Default: the isolates being designed")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--min-hamming", type=int, default=2)
    p.add_argument("--num-barcodes", type=int, default=100,
                   help="barcodes sampled per target in Step 0")
    p.add_argument("--flank-size", type=int, default=300)
    p.add_argument("--product-size-range", default="80,120",
                   help="Primer3 PRIMER_PRODUCT_SIZE_RANGE, e.g. '80,120'")
    p.add_argument("--retry", action="store_true",
                   help="multi-round retry loop instead of a single pass")
    p.add_argument("--rounds", type=int, default=10)
    p.add_argument("--use-existing-candidates", action="store_true",
                   help="also load Stage 2 CSVs (only safe if regenerated against "
                        "the current reference set)")
    p.add_argument("--no-fresh-candidates", action="store_true",
                   help="skip Primer3 and use Stage 2 CSVs only")
    p.add_argument("--skip-locked-audit", action="store_true")
    p.add_argument("--skip-locked-pair-audit", action="store_true",
                   help="skip check A3 only (the slow one)")
    p.add_argument("--max-pcr-product", type=int, default=300)
    p.add_argument("--seed-len", type=int, default=9,
                   help="3' exact seed length for off-target search (Stage 3 uses 9)")
    p.add_argument("--max-mismatches", type=int, default=2)
    return p.parse_args(argv[1:])


def main(argv):
    args = parse_args(argv)
    max_length = args.max_length if args.max_length is not None else args.length
    suffix = (f"k{args.length}_{max_length}" if max_length != args.length
              else f"k{args.length}")

    BARCODE_CSV = f"./output/KMERS/unique_barcodes_{suffix}.csv"
    CANDIDATE_DIR = f"./output/barcodes_{suffix}_amplicon_candidates/"
    OUTPUT_DIR = args.out_dir or f"./output/barcodes_{suffix}_amplicon_addition/"

    np.random.seed(args.seed)
    random.seed(args.seed)

    psr = [int(x) for x in args.product_size_range.split(",")]
    primer_config = {
        "flank_size": args.flank_size,
        "primer3_globals": {
            "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
            "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0,
            "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
            "PRIMER_PRODUCT_SIZE_RANGE": [psr],
        },
    }

    ### Reference Sequences
    print("Loading chromosome templates...")
    target_seqs = parse_sequences(find_fasta_files(args.only_genome))
    if not target_seqs:
        print(f"\nERROR: no FASTA files in {args.only_genome}.")
        return 1

    print("\nLoading background records (specificity reference)...")
    all_records, stem_to_ids = parse_background_by_file(args.background)
    if not all_records:
        print(f"\nERROR: no records in {args.background}.")
        return 1

    print(f"\nLoading locked panel: {args.locked}")
    locked = load_locked_panel(args.locked)
    locked_keys = set(locked.keys())
    print(f"  {len(locked)} locked amplicons "
          f"({len({k[0] for k in locked_keys})} isolates).")

    ### Release amplicons for redesign
    unlocked = []
    if args.unlock:
        for key in load_targets_file(args.unlock):
            if key in locked:
                del locked[key]
                unlocked.append(key)
            else:
                print(f"  NOTE: --unlock entry {_fmt(key)} is not in the locked panel")
        if unlocked:
            print(f"  Released {len(unlocked)} amplicons for re-design: "
                  f"{[_fmt(k) for k in unlocked]}")
            print("  Their previously synthesized primers will NOT appear in the output.")
        locked_keys = set(locked.keys())

    if args.targets:
        targets = load_targets_file(args.targets)
    else:
        targets = [(sid, r) for sid in sorted(target_seqs) for r in REGIONS]
    targets = list(unlocked) + [t for t in targets if t not in locked_keys]
    # de-duplicate, preserve order
    seen = set()
    targets = [t for t in targets if not (t in seen or seen.add(t))]

    if args.new_isolates:
        new_seq_ids = [l.strip() for l in open(args.new_isolates)
                       if l.strip() and not l.startswith("#")]
    else:
        new_seq_ids = sorted({sid for sid, _ in targets})
    new_record_ids = []
    for sid in new_seq_ids:
        stem = isolate_stem_for(sid, stem_to_ids)
        new_record_ids.extend(stem_to_ids.get(stem, []))
    new_record_ids = sorted(set(new_record_ids))

    # Audit
    if not args.skip_locked_audit:
        print(f"\nStep A: auditing the locked panel...")
        t0 = timer.perf_counter()
        report = audit_locked_panel(
            locked, all_records, new_record_ids,
            max_pcr_product=args.max_pcr_product, seed_len=args.seed_len,
            max_mismatches=args.max_mismatches,
            check_locked_pairs=not args.skip_locked_pair_audit,
        )
        print_audit(report, len(locked), len(new_record_ids))
        print(f"  (audit took {timer.perf_counter() - t0:.1f}s)")
        write_audit(report, OUTPUT_DIR, suffix)

    if args.mode == "audit":
        return 0

    print(f"\n{len(targets)} targets to design.")
    if not targets:
        print("Nothing to do — every isolate/region is already in the locked panel.")
        return 0

    # Builds Candidates
    print(f"\nStep 0: building candidate pools...")
    t0 = timer.perf_counter()
    existing = load_all_candidate_primers(CANDIDATE_DIR) if args.use_existing_candidates else {}
    new_candidates = {}
    no_candidates = []
    exhausted_pool = []

    for key in targets:
        seq_id, region = key
        old = existing.get(key, [])
        fresh = []
        n_barcodes = 0

        if not args.no_fresh_candidates:
            if seq_id not in target_seqs:
                print(f"  {_fmt(key)}: no chromosome FASTA in {args.only_genome}; skipped")
                no_candidates.append(key)
                continue
            barcodes = load_unique_barcodes(BARCODE_CSV, seq_id, region=region or None)
            n_barcodes = len(barcodes)
            if barcodes:
                if n_barcodes <= args.num_barcodes:
                    exhausted_pool.append(key)
                fresh = design_candidate_primers(
                    target_seqs, all_records, seq_id, barcodes,
                    primer_config, args.num_barcodes,
                )
            else:
                print(f"  {_fmt(key)}: no unique barcodes in {BARCODE_CSV}")

        combined = old + fresh
        pool = "" if args.no_fresh_candidates else f" (from {n_barcodes} unique barcodes)"
        if combined:
            new_candidates[key] = combined
            flag = " [pool exhausted]" if key in exhausted_pool else ""
            print(f"  {_fmt(key)}: {len(old)} existing + {len(fresh)} fresh "
                  f"= {len(combined)}{pool}{flag}")
        else:
            no_candidates.append(key)
            print(f"  {_fmt(key)}: no candidates{pool}")

    print(f"  Candidate generation done ({timer.perf_counter() - t0:.1f}s)")

    if exhausted_pool:
        print(f"\n  NOTE: {len(exhausted_pool)} targets have <= --num-barcodes "
              f"({args.num_barcodes}) unique barcodes, so every seed samples their "
              f"entire barcode space and builds an identical candidate pool:")
        print(f"    {[_fmt(k) for k in sorted(exhausted_pool)]}")
        print("  More seeds cannot help these — change the barcode length, widen "
              "--product-size-range, or relax the Primer3 limits.")

    if not new_candidates:
        print("\nNo candidates available. Nothing to add.")
        save_panel(OUTPUT_DIR, locked, f"added_primers_{suffix}-seed{args.seed}.csv")
        return 0

    # Pre-filter
    print(f"\nStep 1: pre-filtering against {len(locked)} locked amplicons...")
    t0 = timer.perf_counter()
    filtered, prefilter_failed = prefilter_candidates(
        locked, new_candidates, args.min_hamming)
    before = sum(len(v) for v in new_candidates.values())
    after = sum(len(v) for v in filtered.values())
    print(f"  {before} -> {after} candidates compatible with the locked panel "
          f"({timer.perf_counter() - t0:.1f}s)")
    print(f"  {len(filtered)} targets retained, {len(prefilter_failed)} eliminated")

    if filtered:
        print(f"\nSteps 2-4: selection and cross-reactivity "
              f"({'retry loop' if args.retry else 'single pass'})...")
        accepted, selection_failed, exhausted, still_pending = select_and_validate(
            filtered, locked, all_records, args.min_hamming,
            max_rounds=args.rounds if args.retry else 1,
            max_pcr_product=args.max_pcr_product, seed_len=args.seed_len,
            max_mismatches=args.max_mismatches,
        )
    else:
        accepted, selection_failed, exhausted, still_pending = {}, set(), set(), set()

    merged = {**locked, **accepted}
    still_failed = (set(no_candidates) | set(prefilter_failed) |
                    selection_failed | exhausted | still_pending) - set(accepted)

    print(f"\n--- SUMMARY (seed {args.seed}) ---")
    print(f"  Locked (unchanged):  {len(locked)}")
    print(f"  Newly designed:      {len(accepted)}")
    print(f"  Final panel size:    {len(merged)}")
    print(f"  Unresolved targets:  {len(still_failed)}")

    save_panel(OUTPUT_DIR, merged, f"added_primers_{suffix}-seed{args.seed}.csv")
    if accepted:
        save_panel(OUTPUT_DIR, accepted,
                   f"added_primers_NEW_ONLY_{suffix}-seed{args.seed}.csv")

    if still_failed:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        path = os.path.join(OUTPUT_DIR, f"still_failed_{suffix}-seed{args.seed}.txt")
        with open(path, "w") as f:
            for seq_id, region in sorted(still_failed):
                f.write(f"{seq_id},{region}\n")
        print(f"  Unresolved targets -> {path}")

        print("\n--- FAILURE BREAKDOWN ---")
        if no_candidates:
            print(f"  No candidates at all ({len(no_candidates)}): "
                  f"{[_fmt(k) for k in sorted(no_candidates)]}")
        if prefilter_failed:
            print(f"  Incompatible with locked panel ({len(prefilter_failed)}): "
                  f"{[_fmt(k) for k in sorted(prefilter_failed)]}")
        sel = sorted(selection_failed - exhausted)
        if sel:
            print(f"  Mutually incompatible among new targets ({len(sel)}): "
                  f"{[_fmt(k) for k in sel]}")
        if exhausted:
            print(f"  Cross-reactivity, candidates exhausted ({len(exhausted)}): "
                  f"{[_fmt(k) for k in sorted(exhausted)]}")

    print("\n--- COMPLETE ---")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
