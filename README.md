# multiplex_PCR_design

Authors: Jaemin Kim, Zhengqing Zhou. (Duke University)

Find a unique 20 bp sequence within each strain's Ori and Ter regions, then design and validate PCR primers around each one so all 36 amplicons can be distinguished from each other and amplified together in a single multiplex reaction without cross-reacting.

---

## Directory layout

```
project/
├── barcode_ranges_in_genome.csv          # Ori_5%_Range and Ter_5%_Range per isolate
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
