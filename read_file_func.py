"""
Utility functions for reading files.
"""

import os
import csv
import fnmatch
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


### FASTA discovery and parsing

def find_fasta_files(input_dir, pattern=None):
    """
    Find supported sequence files in a directory.

    If `pattern` is provided (e.g. "*_p*.fasta"), only return matching files.
    If `pattern` is None, return all files with supported extensions.
    Skips macOS metadata files (names starting with "._").
    """
    print("Searching for sequence files...")
    fasta_files = []
    supported_extensions = (".fasta", ".fa", ".fna", ".txt")

    for filename in os.listdir(input_dir):
        if filename.startswith("._"):
            continue
        full_path = os.path.join(input_dir, filename)
        if not os.path.isfile(full_path):
            continue
        if pattern is not None:
            if fnmatch.fnmatch(filename.lower(), pattern.lower()):
                fasta_files.append(full_path)
        else:
            if filename.lower().endswith(supported_extensions):
                fasta_files.append(full_path)

    if not fasta_files:
        if pattern is not None:
            print(f"Warning: No files matching '{pattern}' found in '{input_dir}'.")
        else:
            print(f"Warning: No files with extensions {supported_extensions} found in '{input_dir}'.")
    else:
        print(f"Found {len(fasta_files)} sequence file(s).")

    return fasta_files


def parse_sequences(fasta_files):
    """
    Parse sequences from a list of FASTA files.
    Each file contributes ONE entry keyed by the filename stem
    (e.g. 'KL13_1' from 'KL13_1.fasta').  When a file has multiple records,
    only the first record is stored — use parse_all_records_from_dir() if you
    need every record (e.g. chromosomes + plasmids).
    """
    print("Parsing sequences...")
    sequences = {}
    for file_path in fasta_files:
        try:
            seq_id  = os.path.splitext(os.path.basename(file_path))[0]
            records = list(SeqIO.parse(file_path, "fasta-pearson"))
            if records:
                print(f"  {file_path} — OK ({len(records)} record(s), using first)")
                sequences[seq_id] = records[0].seq
            else:
                print(f"  {file_path} — FASTA parse failed; trying raw text")
                with open(file_path, 'r') as f:
                    raw = "".join(
                        line.strip() for line in f
                        if line.strip() and not line.startswith(">")
                    )
                if raw:
                    sequences[seq_id] = Seq(raw)
        except Exception as e:
            print(f"  Error processing {file_path}: {e}")
    print(f"\nSuccessfully parsed {len(sequences)} sequences.")
    return sequences


def parse_all_records_from_dir(directory: str) -> dict[str, str]:
    """
    Load EVERY FASTA record from EVERY file in `directory`.
    Returns {record.id: sequence_str}.
    Used to build a complete background for specificity checks (chromosomes +
    plasmids for all isolates).

    Skips macOS metadata files (names starting with '._').
    """
    records_dict: dict[str, str] = {}
    dir_path = Path(directory)
    fasta_files = sorted(
        [f for f in dir_path.glob("*.fasta") if not f.name.startswith("._")] +
        [f for f in dir_path.glob("*.fa")    if not f.name.startswith("._")]
    )
    for fasta_file in fasta_files:
        for rec in SeqIO.parse(str(fasta_file), "fasta"):
            records_dict[rec.id] = str(rec.seq).upper()
    print(f"Loaded {len(records_dict)} records from {len(fasta_files)} files in '{directory}'.")
    return records_dict


def parse_all_records_except(directory: str, exclude_prefix: str) -> dict[str, str]:
    """
    Same as parse_all_records_from_dir() but skips files whose stem starts
    with `exclude_prefix` (e.g. exclude the target isolate's own file so we
    don't flag the intended amplification site as off-target).

    Returns {record.id: sequence_str}.
    """
    records_dict: dict[str, str] = {}
    dir_path = Path(directory)
    fasta_files = sorted(
        [f for f in dir_path.glob("*.fasta") if not f.name.startswith("._")] +
        [f for f in dir_path.glob("*.fa")    if not f.name.startswith("._")]
    )
    skipped = 0
    for fasta_file in fasta_files:
        stem = fasta_file.stem  # filename without extension
        if stem.startswith(exclude_prefix):
            skipped += 1
            continue
        for rec in SeqIO.parse(str(fasta_file), "fasta"):
            records_dict[rec.id] = str(rec.seq).upper()
    print(f"Loaded {len(records_dict)} background records "
          f"(skipped {skipped} file(s) matching '{exclude_prefix}*').")
    return records_dict



### Barcode CSV  (output of kmer_generation.py)

def load_unique_barcodes(csv_path: str, seq_id: str,
                         region: str | None = None) -> list[tuple[str, int]]:
    """
    Load barcode candidates for a specific seq_id from the kmer_generation
    output CSV (which may contain barcodes of multiple lengths).

    Returns:
        [(barcode_sequence, position_0based), ...]  sorted by position.
        Barcodes of all lengths within the file are included; the primer
        design stage handles variable-length barcodes correctly.
        Position is 0-based chromosomal coordinate.
    """
    barcodes = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["Seq_ID"] != seq_id:
                continue
            if region is not None and row.get("Region", "") != region:
                continue
            barcodes.append((row["Barcode_Sequence"], int(row["Position"])))
    return sorted(barcodes, key=lambda x: x[1])


### Ranges CSV  (barcode_ranges_in_genome.csv)

def parse_ranges_csv(csv_path: str) -> list[dict]:
    """
    Parse barcode_ranges_in_genome.csv.
    Returns a list of dicts, one per isolate, with column names as keys.
    Column names are stripped of leading/trailing whitespace.
    """
    rows = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        # strip any whitespace from column headers
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})
    return rows


### Candidate primer CSVs  (output of individual_target_amplicon_candidates.py)

def load_candidate_primers(csv_path: str) -> dict:
    """
    Load a single candidate primer CSV.
    Filename format: SEQ_ID_Ori.csv or SEQ_ID_Ter.csv
    Returns {(seq_id, region): [candidate_dict, ...]}.
    """
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    # Parse region from end of stem: e.g. "KL13_1_Ori" → seq_id="KL13_1", region="Ori"
    if stem.endswith('_Ori'):
        seq_id = stem[:-4]
        region = 'Ori'
    elif stem.endswith('_Ter'):
        seq_id = stem[:-4]
        region = 'Ter'
    else:
        # Fallback for old-style files without region suffix
        seq_id = stem
        region = 'Unknown'

    candidates = []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append({
                "target_barcode":   row["Target_Barcode"],
                "fwd_primer":       row["Forward_Primer"],
                "rev_primer":       row["Reverse_Primer"],
                "product_size":     int(row["Product_Size"]),
                "fwd_tm":           float(row["Fwd_Tm"]),
                "rev_tm":           float(row["Rev_Tm"]),
                "amplicon_sequence": row["Amplicon_Sequence"],
            })
    return {(seq_id, region): candidates}


def load_all_candidate_primers(csv_dir: str) -> dict:
    """
    Load every *_Ori.csv and *_Ter.csv in csv_dir.
    Returns {(seq_id, region): [candidate_dict, ...]}.
    """
    all_candidates = {}
    for fname in os.listdir(csv_dir):
        if not fname.endswith(".csv"):
            continue
        all_candidates.update(load_candidate_primers(os.path.join(csv_dir, fname)))
    return all_candidates
