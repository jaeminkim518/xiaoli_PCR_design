"""
Find candidate amplicons and primers for each sequence.

Usage (called from SLURM array script):
    python individual_target_amplicon_candidates.py \\
        <RANGES_CSV> <ONLY_GENOME_DIR> <BACKGROUND_DIR> <SEQ_ID> <BARCODE_LENGTH>
"""

import sys
import os
import re
# Ensure the script's own directory is on the path so local modules are found
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import primer3
import numpy as np
import time as timer

from read_file_func import (
    find_fasta_files,
    parse_sequences,
    parse_all_records_from_dir,
    parse_all_records_except,
    load_unique_barcodes,
)
from sequence_alignment import is_primer_pair_specific


### Primer Design
def design_candidate_primers(
    target_seqs: dict,           # {SEQ_ID: chrom_seq}   — template source
    specificity_seqs: dict,      # {any_id: seq_str}      — off-target check
    seq_id: str,
    barcodes: list[tuple[str, int]],
    primer_config: dict,
    num_of_barcodes_checked: int = 50,
) -> list[dict]:
    """
    Design primer pairs that flank each candidate barcode within the target
    chromosome and are verified not to amplify anywhere in specificity_seqs.
    """
    candidates    = []
    target_seq    = str(target_seqs[seq_id])

    # Random subsample if more barcodes than we can afford to check
    if len(barcodes) <= num_of_barcodes_checked:
        candidate_barcodes = barcodes
    else:
        idxs = np.random.choice(len(barcodes), num_of_barcodes_checked, replace=False)
        candidate_barcodes = [barcodes[i] for i in idxs]

    for barcode, position in candidate_barcodes:
        flank_size    = primer_config['flank_size']
        template_start = max(0, position - flank_size)
        template_end   = min(len(target_seq), position + len(barcode) + flank_size)
        sequence_template        = target_seq[template_start:template_end]
        target_start_in_template = position - template_start

        primer3_results = primer3.design_primers(
            seq_args={
                'SEQUENCE_TEMPLATE': sequence_template,
                'SEQUENCE_TARGET':   [target_start_in_template, len(barcode)],
            },
            global_args=primer_config['primer3_globals'],
        )

        for i in range(primer3_results['PRIMER_PAIR_NUM_RETURNED']):
            fwd_primer = primer3_results[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev_primer = primer3_results[f'PRIMER_RIGHT_{i}_SEQUENCE']

            # Specificity check against target + all off-target sequences
            if is_primer_pair_specific(fwd_primer, rev_primer, seq_id, specificity_seqs):
                left_start,  _ = primer3_results[f'PRIMER_LEFT_{i}']
                right_start, _ = primer3_results[f'PRIMER_RIGHT_{i}']
                amplicon_seq   = sequence_template[left_start:right_start + 1]
                prod           = primer3_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                assert (right_start - left_start + 1) == prod

                candidates.append({
                    'target_barcode':   barcode,
                    'fwd_primer':       fwd_primer,
                    'rev_primer':       rev_primer,
                    'product_size':     prod,
                    'fwd_tm':           round(primer3_results[f'PRIMER_LEFT_{i}_TM'],  2),
                    'rev_tm':           round(primer3_results[f'PRIMER_RIGHT_{i}_TM'], 2),
                    'amplicon_sequence': amplicon_seq,
                })

    return candidates


def save_candidate_primers(seq_id: str, candidates: list, output_dir: str, filename: str):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    with open(output_path, 'w', newline='') as f:
        header = ["Seq_ID", "Target_Barcode", "Forward_Primer", "Reverse_Primer",
                  "Product_Size", "Fwd_Tm", "Rev_Tm", "Amplicon_Sequence"]
        f.write(",".join(header) + "\n")
        for data in candidates:
            row = [
                str(seq_id),
                data['target_barcode'],
                data['fwd_primer'],
                data['rev_primer'],
                str(data['product_size']),
                str(data['fwd_tm']),
                str(data['rev_tm']),
                data['amplicon_sequence'],
            ]
            f.write(",".join(row) + "\n")
    print(f"Saved {len(candidates)} candidates for {seq_id} → {output_path}")


if __name__ == "__main__":
    np.random.seed(42)

    if len(sys.argv) < 6:
        print("Usage: python individual_target_amplicon_candidates.py "
              "<RANGES_CSV> <ONLY_GENOME_DIR> <BACKGROUND_DIR> <SEQ_ID> <LENGTH> [MAX_LENGTH]")
        sys.exit(1)

    RANGES_CSV      = sys.argv[1]
    ONLY_GENOME_DIR = sys.argv[2]
    BACKGROUND_DIR  = sys.argv[3]
    SEQ_ID          = sys.argv[4]
    MIN_LENGTH      = int(sys.argv[5])
    MAX_LENGTH      = int(sys.argv[6]) if len(sys.argv) > 6 else MIN_LENGTH

    barcode_suffix = f"k{MIN_LENGTH}_{MAX_LENGTH}" if MAX_LENGTH != MIN_LENGTH else f"k{MIN_LENGTH}"
    BARCODE_CSV = f"./output/KMERS/unique_barcodes_{barcode_suffix}.csv"
    OUTPUT_DIR  = f"./output/barcodes_{barcode_suffix}_amplicon_candidates/"

    # Target chromosome
    target_seqs_all = parse_sequences(find_fasta_files(ONLY_GENOME_DIR))
    if SEQ_ID not in target_seqs_all:
        print(f"\nERROR: '{SEQ_ID}' not found in {ONLY_GENOME_DIR}.")
        sys.exit(1)
    target_seqs = {SEQ_ID: target_seqs_all[SEQ_ID]}

    # Off-target background
    off_target_seqs = parse_all_records_from_dir(BACKGROUND_DIR)
    specificity_seqs = off_target_seqs

    # Primer3 config
    PRODUCT_SIZE_RANGE = [80, 120]
    primer_config = {
        'flank_size': 300,
        'primer3_globals': {
            'PRIMER_OPT_SIZE':           20,
            'PRIMER_MIN_SIZE':           18,
            'PRIMER_MAX_SIZE':           25,
            'PRIMER_OPT_TM':            60.0,
            'PRIMER_MIN_TM':            57.0,
            'PRIMER_MAX_TM':            63.0,
            'PRIMER_MIN_GC':            40.0,
            'PRIMER_MAX_GC':            60.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [PRODUCT_SIZE_RANGE],
        },
    }

    t0 = timer.perf_counter()
    num_of_barcodes_checked = 100

    # Design candidates separately for Ori and Ter
    for region in ['Ori', 'Ter']:
        print(f"\n{'='*55}")
        print(f"Region: {region} — {SEQ_ID}")

        barcodes = load_unique_barcodes(BARCODE_CSV, SEQ_ID, region=region)
        if not barcodes:
            print(f"  No barcodes found for {SEQ_ID} {region} in {BARCODE_CSV}.")
            continue
        print(f"  Loaded {len(barcodes)} barcode candidates.")

        candidates = design_candidate_primers(
            target_seqs, specificity_seqs,
            SEQ_ID, barcodes, primer_config,
            num_of_barcodes_checked,
        )
        if not candidates:
            print(f"  Could not design any primers for {SEQ_ID} {region}.")
            continue

        filename = f"{SEQ_ID}_{region}.csv"
        save_candidate_primers(SEQ_ID, candidates, OUTPUT_DIR, filename=filename)

    print(f"\n{(timer.perf_counter() - t0) / 60:.1f} min elapsed.")
