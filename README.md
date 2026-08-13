# multiplex_PCR_design

Authors: Jaemin Kim, Zhengqing Zhou. (Duke University)

Find a unique 20 bp sequence within each strain's Ori and Ter regions, then design and validate PCR primers around each one so all 36 amplicons can be distinguished from each other and amplified together in a single multiplex reaction without cross-reacting.

---

## Directory layout

```
project/
├── barcode_ranges_in_genome.csv          # Ori_5%_Range and Ter_5%_Range per isolate
├── amplicon_design_env.yml
├── only_genome/                          # chromosome-only FASTA, one file per isolate
│   ├── KL1_1.fasta
│   ├── KL13_1.fasta
│   └── ...                              # naming convention: <isolate>_1.fasta
├── genome_and_plasmids_within_host/      # chromosome + plasmids, one file per isolate
│   ├── KL1.fasta
│   ├── KL13.fasta
│   └── ...                              # naming convention: <isolate>.fasta
├── kmer_generation.py
├── individual_target_amplicon_candidates.py
├── design_amplicon.py
├── read_file_func.py
├── sequence_alignment.py
├── candidates.sh
└── design_array.sh
└── add_amplicons.py
└── addition_array.sh                 # writes to ./output/barcodes_k20_amplicon_addition/
```

---

## Installation (one-time)

```bash
conda create -n amplicon_design -f amplicon_design_env.yml
conda activate amplicon_design
```

---

## Running the pipeline

The pipeline has three stages that must be run in order.
Stages 2 and 3 are parallelized via SLURM array jobs.

Stage 4 is optional. Run it after a part of the panel is already fixed - if new isolates need adding, or if a Stage 3 run left some targets without an amplicon.

### Stage 1 — Find unique barcodes *(run directly on compute node, ~5–30 min)*

Get an interactive compute node:
```bash
srun -p youlab-gpu --cpus-per-task=2 --mem=128G -t 1-00:00:00 --pty /bin/bash
```

Activate the environment and run:
```bash
conda activate amplicon_design
cd ~/xiaoli_PCR_design

python kmer_generation.py \
    barcode_ranges_in_genome.csv \
    ./only_genome \
    ./genome_and_plasmids_within_host \
    20
```

Output: `./output/KMERS/unique_barcodes_k20.csv`

This file contains all unique 20 bp barcode candidates for every strain's
Ori and Ter regions, verified unique against all 18 genomes and all plasmids.

---

### Stage 2 — Design primer candidates *(SLURM array, one job per isolate)*

Before submitting, check the number of isolates and confirm the `--array`
bound in `candidates.sh` matches (0 to N-1):
```bash
ls ./only_genome/*.fasta | grep -v '\._' | wc -l
```

Submit:
```bash
sbatch candidates.sh
```

Monitor until all jobs complete:
```bash
squeue -u <your_username>
```

Output: `./output/barcodes_k20_amplicon_candidates/`
One file per strain per region: `SEQ_ID_Ori.csv` and `SEQ_ID_Ter.csv`

Each file contains primer pairs that:
- Flank a unique barcode within the correct Ori or Ter region
- Pass Primer3 checks (GC 40-60%, Tm 57-63°C, length 18-25 bp)
- Do not amplify any other strain's genome or plasmids

Do NOT proceed to Stage 3 until all Stage 2 jobs are finished.

---

### Stage 3 — Select optimal multiplex combination *(SLURM array, one job per seed)*

```bash
sbatch design_array.sh
```

Monitor until all jobs complete:
```bash
squeue -u <your_username>
```

Output: `./output/barcodes_k20_amplicon_design/`
One CSV per seed: `validated_primers_k20-seed<N>.csv`

Each output file contains up to 36 primer pairs (18 Ori + 18 Ter) that:
- Are collectively unique and distinguishable (Hamming distance check)
- Do not form primer dimers with each other
- Do not cross-amplify any other strain when run together in multiplex

---

### Selecting the best result

Find the seed with the most amplicons:
```bash
wc -l ./output/barcodes_k20_amplicon_design/*.csv
```
The file with the most lines (subtract 1 for the header) is your best result.
A perfect run yields 37 lines (36 amplicons + 1 header).

---

### Stage 3 — Select optimal multiplex combination *(SLURM array, one job per seed)*

```bash
sbatch design_array.sh
```

Monitor until all jobs complete:
```bash
squeue -u <your_username>
```

Output: `./output/barcodes_k20_amplicon_design/`
One CSV per seed: `validated_primers_k20-seed<N>.csv`

Each output file contains up to 36 primer pairs (18 Ori + 18 Ter) that:
- Are collectively unique and distinguishable (Hamming distance check)
- Do not form primer dimers with each other
- Do not cross-amplify any other strain when run together in multiplex

---

### Selecting the best result

Find the seed with the most amplicons:
```bash
wc -l ./output/barcodes_k20_amplicon_design/*.csv
```
The file with the most lines (subtract 1 for the header) is your best result.
A perfect run yields 37 lines (36 amplicons + 1 header).

---

### Stage 4 — Add amplicons to an existing panel *(optional)*

Everything in the locked CSV is copied to the output verbatim; only the targets
you ask for are designed, and they are designed to be compatible with what is
already there. Targets are `(Seq_ID, Region)` pairs, the same keys Stage 3 uses.

#### Adding new isolates

**1. Add the new isolates to the existing directories.** No new directories.

```bash
cp KL30_1.fasta ./only_genome/                        # chromosome only
cp KL30.fasta   ./genome_and_plasmids_within_host/    # chromosome + plasmids
```

The chromosome appears in both directories, which is intended: `only_genome/`
supplies the primer-design template, `genome_and_plasmids_within_host/` is the
off-target reference. The `only_genome` filename stem (`KL30_1`) becomes the
`Seq_ID` used downstream and must match the chromosome's FASTA record id.

Then add a row for the new isolate to `barcode_ranges_in_genome.csv`, giving the
5% windows around its origin and terminus as 1-based inclusive coordinates:

```
isolate,Ori_5%_Range,Ter_5%_Range
KL30,1 - 250000,2400000 - 2650000
```

The `isolate` column is `KL30`, not `KL30_1` — Stage 1 maps it to the file by
prefix. A window that wraps the origin is written as two intervals joined by
`and`, e.g. `6666100 - 7010776 and 0 - 356401`. Without this row Stage 1 skips
the isolate and finds no barcodes for it.

You do **not** need to re-run Stage 2 or touch the `--array` bound in
`candidates.sh`: Stage 4 designs its own candidates for the new targets. Only do
so if you want the Stage 2 CSVs regenerated for some other reason.

**2. Re-run Stage 1 over the combined set.** Mandatory: a
barcode is only unique relative to the FASTAs present when Stage 1 ran, so new
isolates can retroactively break barcodes that are already synthesized.

```bash
python kmer_generation.py \
    barcode_ranges_in_genome.csv \
    ./only_genome \
    ./genome_and_plasmids_within_host \
    20
```

**3. Audit the locked panel** before designing anything:

```bash
python add_amplicons.py --mode audit \
    --locked ./output/barcodes_k20_amplicon_design/validated_primers_k20-seed7.csv \
    --length 20
```

**4. Design the new targets.** With no `--targets`, every isolate/region pair in
`only_genome/` that is not already in the locked CSV is designed:

```bash
python add_amplicons.py \
    --locked ./output/barcodes_k20_amplicon_design/validated_primers_k20-seed7.csv \
    --length 20 --seed 1
```

Across seeds (the locked rows are identical in every output; only the new
amplicons differ):

```bash
sbatch addition_array.sh          # edit LOCKED_CSV at the top first
wc -l ./output/barcodes_k20_amplicon_addition/added_primers_k20-seed*.csv
```

#### Re-designing specific targets

To retry targets a Stage 3 run could not place, name them explicitly. Either
form works, one per line: `SEQ_ID,Region`, or a bare `SEQ_ID` for both regions.

```bash
printf 'KL13_1,Ori\nKL13_1,Ter\n' > targets.txt
python add_amplicons.py --locked <csv> --length 20 --seed 1 \
    --targets targets.txt --skip-locked-audit
```

To replace an amplicon that is *already in* the locked panel — for example one
the audit flagged — release it with `--unlock`. Its old primers will not appear
in the output, so only do this for amplicons you are prepared to re-order.

```bash
printf 'KL1_1,Ter\n' > unlock.txt
python add_amplicons.py --locked <csv> --length 20 --seed 1 --unlock unlock.txt
```

#### What the audit checks

Stages 1–3 guarantee uniqueness and specificity only against the reference set
as it stood when they ran. Adding isolates can break three things, none of which
Stages 1–3 re-check:

- **A1** — a locked barcode now also occurs in a new isolate's genome or plasmid
- **A2** — a locked primer pair now amplifies a new record off-target
- **A3** — two locked primer pairs cross-amplify on a new record (forward from
  one, reverse from another, within 300 bp)

A preliminary check (**A0**) confirms each locked `Seq_ID` matches a record id in
the background directory; if it does not, the pair's own amplicon cannot be
excluded and A1/A2 results for it are meaningless rather than merely wrong.

The audit **reports and never repairs** — the oligos already exist. Findings go
to `locked_audit_k20.txt`; an empty file means the panel is clean. A3 is the slow
check and is restricted to the newly added records; `--skip-locked-pair-audit`
turns it off, `--skip-locked-audit` skips Step A entirely (safe when nothing was
added to the reference set).

#### Outputs

Written to `./output/barcodes_k20_amplicon_addition/`:

| File | Contents |
|---|---|
| `locked_audit_k20.txt` | Audit findings, tab-separated; empty means clean |
| `added_primers_k20-seed<N>.csv` | Full panel: locked rows verbatim + new rows, same 9 columns as `validated_primers_*.csv` |
| `added_primers_NEW_ONLY_k20-seed<N>.csv` | Only the new rows — the list to send for synthesis |
| `still_failed_k20-seed<N>.txt` | `SEQ_ID,Region` per line for targets with no compatible amplicon |

Confirm the locked rows really are untouched before ordering:

```bash
tail -n +2 ./output/barcodes_k20_amplicon_design/validated_primers_k20-seed7.csv | sort > /tmp/old.txt
grep -Ff <(cut -d, -f1 ./output/barcodes_k20_amplicon_design/validated_primers_k20-seed7.csv | tail -n +2 | sort -u) \
    ./output/barcodes_k20_amplicon_addition/added_primers_k20-seed1.csv | sort > /tmp/new.txt
diff /tmp/old.txt /tmp/new.txt && echo "locked panel unchanged"
```

#### Notes

- Step 0 designs fresh candidates with Primer3 and re-checks specificity against
  the current reference, so **do not** reuse Stage 2 CSVs written before the new
  isolates were added — they were validated against the smaller set.
  `--use-existing-candidates` is opt-in for exactly that reason.
- Each seed samples `--num-barcodes` (default 100) barcodes per target, so
  different seeds explore different candidate pools. When a target has fewer
  unique barcodes than that, every seed uses all of them and the pool is
  identical; Step 0 flags these as `[pool exhausted]`, meaning more seeds cannot
  help and the barcode length or size range has to change instead.
- `--retry` re-runs selection for up to 10 rounds, blacklisting the candidate
  that failed cross-reactivity. It costs roughly 6x the runtime and rarely
  changes the outcome — prefer more seeds with the default single pass.

---

## Key parameters

All parameters are set at the top of each script or shell file.

| Parameter | Location | Default | Description |
|---|---|---|---|
| Barcode length | `candidates.sh`, `design_array.sh` | 20 bp | Set MIN_LENGTH and MAX_LENGTH (equal for single length) |
| PRODUCT_SIZE_RANGE | `individual_target_amplicon_candidates.py` | [80, 120] | Total amplicon size in bp (barcode + both primers) |
| num_of_barcodes_checked | `individual_target_amplicon_candidates.py` | 100 | Barcodes sampled per strain per region in Stage 2 |
| min_hamming_dist | `design_amplicon.py` | 2 | Min edit distance between same-length barcodes in Stage 3 |
| --array | `candidates.sh` | 0-17 | Must match number of isolates (0 to N-1) |
| --array | `design_array.sh` | 1-20 | Number of random seeds to try |
