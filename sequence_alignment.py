from Bio.Seq import Seq
from bisect import bisect_left

def window_passes(window: str, primer: str, max_mismatches: int) -> bool:
    mism = 0
    for a, b in zip(window, primer):
        if a != b:
            mism += 1
            if mism > max_mismatches:
                return False
    return True

def _find_sites_seeded(seq: str, primer: str, seed_len: int, max_mismatches: int):
    """
    Find binding start positions for `primer` in `seq` using:
      1) exact 3' seed match (last `seed_len` bases of primer)
      2) full-length mismatch verification with <= max_mismatches
    """
    L = len(primer)
    if L == 0 or len(seq) < L:
        return []

    seed_len = min(seed_len, L)
    seed = primer[-seed_len:]     # 3' seed (exact)
    offset = L - seed_len         # seed_start = primer_start + offset

    sites = []
    pos = seq.find(seed)
    while pos != -1:
        start = pos - offset      # candidate full-primer start in seq
        if 0 <= start <= len(seq) - L:
            window = seq[start:start + L]
            if window_passes(window, primer, max_mismatches):
                sites.append(start)
        pos = seq.find(seed, pos + 1)

    return sites

def _has_amplicon_on_template(
    seq_str: str,
    left_primer: str,
    right_primer_rc: str,
    max_pcr_product: int,
    seed_len: int,
    max_mismatches: int,
) -> bool:
    """
    Returns True if there exists an oriented amplicon on `seq_str` where:
      left_primer binds at f
      right_primer_rc binds downstream at r
      product_size = (r + len(right_primer_rc)) - f is in (0, max_pcr_product)
    """
    right_len = len(right_primer_rc)

    left_sites = _find_sites_seeded(seq_str, left_primer, seed_len, max_mismatches)
    if not left_sites:
        return False

    right_sites = _find_sites_seeded(seq_str, right_primer_rc, seed_len, max_mismatches)
    if not right_sites:
        return False

    left_sites.sort()
    right_sites.sort()

    for f in left_sites:
        j = bisect_left(right_sites, f + 1)  #   # find the first rev binding site that's downstream of f
        if j == len(right_sites):
            continue
        r = right_sites[j]
        product_size = (r + right_len) - f
        if 0 < product_size < max_pcr_product:
            return True

    return False

def is_primer_pair_specific(
    fwd_primer: str,
    rev_primer: str,
    target_id: str,
    all_sequences: dict,
    max_pcr_product: int = 300,
    seed_len: int = 9,
    max_mismatches: int = 2,
) -> bool:
    """
    Returns True if no off-target sequence contains a plausible amplicon
    in EITHER orientation (i.e., considering both strands).

    We test two orientations on the stored sequence string:
      1) fwd + rev_rc
      2) rev + fwd_rc
    """
    fwd = fwd_primer.upper()
    rev = rev_primer.upper()
    fwd_rc = str(Seq(fwd_primer).reverse_complement()).upper()
    rev_rc = str(Seq(rev_primer).reverse_complement()).upper()

    for seq_id, seq in all_sequences.items():
        if seq_id == target_id:
            continue

        seq_str = str(seq).upper()

        # Orientation 1: fwd binds, rev binds on opposite strand -> search for rev_rc
        if _has_amplicon_on_template(
            seq_str, fwd, rev_rc, max_pcr_product, seed_len, max_mismatches
        ):
            return False

        # Orientation 2: rev binds, fwd binds on opposite strand -> search for fwd_rc
        if _has_amplicon_on_template(
            seq_str, rev, fwd_rc, max_pcr_product, seed_len, max_mismatches
        ):
            return False

    return True