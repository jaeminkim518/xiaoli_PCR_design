"""
Find unique barcode k-mers within the Ori and Ter regions of each strain's
chromosome, verified unique against every genome and plasmid in the dataset.

Usage:
    python kmer_generation.py <RANGES_CSV> <ONLY_GENOME_DIR> \\
                               <BACKGROUND_DIR> <MIN_LENGTH> <MAX_LENGTH>

Output:
    ./output/KMERS/unique_barcodes_k<MIN>_<MAX>.csv  — all lengths combined
    Columns: Seq_ID, Region, Segment, Barcode_Sequence,
             Position, GC_Percent, Barcode_Length
"""

import sys
import os
import re
import time as timer
from pathlib import Path
from collections import defaultdict, Counter

from Bio import SeqIO
from Bio.Seq import Seq

from read_file_func import parse_ranges_csv, find_fasta_files, parse_sequences

### Helper functions

def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _canonical(kmer: str) -> str:
    """Lexicographically smaller of kmer and its reverse complement."""
    rc = _rc(kmer)
    return kmer if kmer <= rc else rc


def parse_range_str(range_str: str) -> list[tuple[int, int]]:
    """
    Parse a range string into a list of (start_1based, end_1based) tuples.
    Handles single intervals and 'and'-joined wrapped intervals.
    Returns each interval separately.
    """
    parts = re.split(r'\band\b', range_str, flags=re.IGNORECASE)
    intervals = []
    for part in parts:
        m = re.match(r'(\d+)\s*-\s*(\d+)', part.strip())
        if m:
            intervals.append((int(m.group(1)), int(m.group(2))))
        else:
            print(f"  [WARN] Could not parse interval: '{part.strip()}'")
    return intervals


def extract_segments(
    chrom_seq: str,
    intervals: list[tuple[int, int]],
) -> list[tuple[str, int, int]]:
    """
    Return one (segment_seq, genomic_start_1based, genomic_end_1based) tuple
    per interval.  CSV coordinates are 1-based inclusive; converted to
    0-based half-open for Python slicing.
    """
    segments = []
    for start_1, end_1 in intervals:
        s = max(0, start_1 - 1)
        e = min(len(chrom_seq), end_1)
        segments.append((chrom_seq[s:e].upper(), start_1, end_1))
    return segments



### Core logic

def find_region_unique_barcodes(
    ranges_rows: list[dict],
    genome_seqs: dict[str, str],        # {seq_id (from filename): chrom_seq}
    isolate_to_seq_id: dict[str, str],  # {csv_isolate_name: seq_id}
    all_background_seqs: list[str],     # every sequence string from background dir
    kmer_size: int,
) -> dict:
    """
    Returns:
        {seq_id: {'Ori': [entry, ...], 'Ter': [entry, ...]}}
        entry keys: region, segment, barcode, position (0-based), gc_percent
    """

    # Slide over Ori/Ter segments and collect candidates
    print(f"\n[Phase 1] Collecting candidates from Ori/Ter segments (k={kmer_size})...")
    candidate_info: dict[str, list] = defaultdict(list)

    for row in ranges_rows:
        isolate = row['isolate']
        seq_id  = isolate_to_seq_id.get(isolate)
        if seq_id is None or seq_id not in genome_seqs:
            print(f"  [SKIP] No genome sequence for isolate '{isolate}'")
            continue
        chrom_seq = genome_seqs[seq_id]

        for region_name, col in [('Ori', 'Ori_5%_Range'), ('Ter', 'Ter_5%_Range')]:
            range_str = str(row.get(col, '')).strip()
            if not range_str or range_str.lower() == 'nan':
                continue

            intervals = parse_range_str(range_str)
            segments  = extract_segments(chrom_seq, intervals)

            for seg_seq, seg_s1, seg_e1 in segments:
                seg_label = f"{seg_s1}-{seg_e1}"
                n = len(seg_seq)
                for i in range(n - kmer_size + 1):
                    kmer = seg_seq[i:i + kmer_size]
                    if 'N' in kmer:
                        continue
                    can   = _canonical(kmer)
                    pos_0 = (seg_s1 - 1) + i   # 0-based chromosomal position
                    candidate_info[can].append(
                        (seq_id, region_name, seg_label, pos_0, kmer)
                    )

    if not candidate_info:
        print("  No candidate k-mers found. Check ranges CSV and genome files.")
        return {}

    candidate_set = set(candidate_info.keys())
    print(f"  {len(candidate_set):,} unique canonical candidates collected.")

    # Count occurrences of candidates in ALL background seqs
    print(f"[Phase 2] Counting across {len(all_background_seqs)} background records...")
    occurrence: Counter = Counter()

    for bg_seq in all_background_seqs:
        n = len(bg_seq)
        for i in range(n - kmer_size + 1):
            kmer = bg_seq[i:i + kmer_size]
            if 'N' in kmer:
                continue
            can = _canonical(kmer)
            if can in candidate_set:
                occurrence[can] += 1

    # Keep only count == 1 (unique across the entire dataset)
    unique_barcodes: dict = defaultdict(lambda: defaultdict(list))

    for can, cnt in occurrence.items():
        if cnt != 1:
            continue
        for (seq_id, region, seg_label, pos_0, original_kmer) in candidate_info[can]:
            gc = (original_kmer.count('G') + original_kmer.count('C')) / kmer_size * 100
            unique_barcodes[seq_id][region].append({
                'region':    region,
                'segment':   seg_label,
                'barcode':   original_kmer,
                'position':  pos_0,        # 0-based, used directly in primer design
                'gc_percent': round(gc, 1),
            })

    # Sort each list by chromosomal position
    for seq_id in unique_barcodes:
        for region in unique_barcodes[seq_id]:
            unique_barcodes[seq_id][region].sort(key=lambda x: x['position'])

    return unique_barcodes


def append_unique_barcodes(
    unique_barcodes: dict,
    kmer_size: int,
    f,          # open file handle, already has header written
) -> int:
    """Append rows for one barcode length to an already-open CSV file."""
    total = 0
    for seq_id in sorted(unique_barcodes.keys()):
        for region in ['Ori', 'Ter']:
            for b in unique_barcodes[seq_id].get(region, []):
                row = [seq_id, b['region'], b['segment'], b['barcode'],
                       str(b['position']), str(b['gc_percent']), str(kmer_size)]
                f.write(",".join(row) + "\n")
                total += 1
    return total


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python kmer_generation.py <RANGES_CSV> <ONLY_GENOME_DIR> "
              "<BACKGROUND_DIR> <LENGTH> [MAX_LENGTH]")
        sys.exit(1)

    RANGES_CSV      = sys.argv[1]
    ONLY_GENOME_DIR = sys.argv[2]
    BACKGROUND_DIR  = sys.argv[3]
    MIN_LENGTH      = int(sys.argv[4])
    MAX_LENGTH      = int(sys.argv[5]) if len(sys.argv) > 5 else MIN_LENGTH

    OUTPUT_DIR      = "./output/KMERS"
    output_filename = (f"unique_barcodes_k{MIN_LENGTH}_{MAX_LENGTH}.csv"
                       if MAX_LENGTH != MIN_LENGTH
                       else f"unique_barcodes_k{MIN_LENGTH}.csv")

    t0 = timer.perf_counter()

    # Ranges CSV
    ranges_rows = parse_ranges_csv(RANGES_CSV)
    print(f"Loaded {len(ranges_rows)} isolates from {RANGES_CSV}.")

    # Chromosomal sequences (only_genome)
    genome_seqs_raw = parse_sequences(find_fasta_files(ONLY_GENOME_DIR))
    genome_seqs     = {sid: str(seq).upper() for sid, seq in genome_seqs_raw.items()}

    # Map CSV isolate name → seq_id (filename stem)
    isolate_to_seq_id: dict[str, str] = {}
    for row in ranges_rows:
        isolate = row['isolate']
        if isolate in genome_seqs:
            isolate_to_seq_id[isolate] = isolate
        else:
            matches = sorted(k for k in genome_seqs if k == isolate or k.startswith(isolate + "_"))
            if matches:
                isolate_to_seq_id[isolate] = matches[0]
                print(f"  Mapped '{isolate}' → '{matches[0]}'")
            else:
                print(f"  [WARN] No genome file found for isolate '{isolate}'")

    # Background: ALL records from genome_and_plasmids dir
    print(f"\nLoading background sequences from: {BACKGROUND_DIR}")
    bg_dir   = Path(BACKGROUND_DIR)
    bg_files = sorted(
        [f for f in bg_dir.glob("*.fasta") if not f.name.startswith("._")] +
        [f for f in bg_dir.glob("*.fa")    if not f.name.startswith("._")]
    )
    all_background_seqs: list[str] = []
    for bg_file in bg_files:
        for rec in SeqIO.parse(str(bg_file), "fasta"):
            all_background_seqs.append(str(rec.seq).upper())
    print(f"  {len(all_background_seqs)} records from {len(bg_files)} files.")

    # Loop over every barcode length in [MIN_LENGTH, MAX_LENGTH]
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path  = os.path.join(OUTPUT_DIR, output_filename)
    grand_total  = 0

    with open(output_path, 'w', newline='') as out_f:
        header = ["Seq_ID", "Region", "Segment", "Barcode_Sequence",
                  "Position", "GC_Percent", "Barcode_Length"]
        out_f.write(",".join(header) + "\n")

        for kmer_size in range(MIN_LENGTH, MAX_LENGTH + 1):
            print(f"\n{'='*55}")
            print(f"  Barcode length: {kmer_size} bp")
            unique_barcodes = find_region_unique_barcodes(
                ranges_rows, genome_seqs, isolate_to_seq_id,
                all_background_seqs, kmer_size,
            )

            if not unique_barcodes:
                print(f"  No unique barcodes found at k={kmer_size}.")
                continue

            n = append_unique_barcodes(unique_barcodes, kmer_size, out_f)
            grand_total += n
            print(f"  k={kmer_size}: {n} barcodes written.")

            # Per-length summary
            for seq_id in sorted(unique_barcodes.keys()):
                ori_n = len(unique_barcodes[seq_id].get('Ori', []))
                ter_n = len(unique_barcodes[seq_id].get('Ter', []))
                print(f"    {seq_id}: Ori={ori_n}, Ter={ter_n}")

    print(f"\n{'='*55}")
    print(f"Done. {grand_total} total barcodes (lengths {MIN_LENGTH}–{MAX_LENGTH} bp)")
    print(f"Output → {output_path}")
    print(f"{(timer.perf_counter() - t0) / 60:.1f} min elapsed.")
