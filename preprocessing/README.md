# Preprocessing

This module converts raw paired-end SELEX FASTQ reads into analysis-ready fixed-length sequence sets.

## What it does

For each sample (two sequencing lanes, paired-end):

1. **Adapter trim + read merge (`fastp`)** — quality filter, TruSeq (or other) adapter removal, and overlap-based merging of R1/R2 into a single read per fragment, independently for lane 1 and lane 2.
2. **Lane concatenation** — concatenate the two lane-merged gzipped FASTQ files into one sample-level FASTQ.gz.
3. **Decompress** — gunzip the combined FASTQ for downstream conversion.
4. **FASTQ → FASTA (`seqtk`)** — convert the combined FASTQ to FASTA.
5. **Length normalization (`filter_by_length.py`)** — keep reads within a user-specified input-length window, then pad or trim so every retained sequence is exactly `--target_length` bases (or the length implied by `--ligand`).
6. **k-mer trim** — extract fixed-length subsequences (default **40-mers**) into a one-sequence-per-line `.seq` file.

## Ligand design expected lengths

| Ligand design | Expected length |
|---------------|-----------------|
| lig147        | 101 bp          |
| lig200        | 154 bp          |
| ligN40        | 187 bp          |
| ligN70        | 247 bp          |

Use `--ligand lig147` (etc.) as a shorthand for the corresponding `--target-length`, or pass an explicit length when needed.

## Generic naming convention

Replace experiment-specific IDs with a sample prefix `SAMPLE`. Expected inputs under `--input-dir`:

| Role | Filename |
|------|----------|
| Lane 1, read 1 | `{SAMPLE}_L1_R1.fastq.gz` |
| Lane 1, read 2 | `{SAMPLE}_L1_R2.fastq.gz` |
| Lane 2, read 1 | `{SAMPLE}_L2_R1.fastq.gz` |
| Lane 2, read 2 | `{SAMPLE}_L2_R2.fastq.gz` |

Outputs written under `--output-dir/{SAMPLE}/` (example with lig147 / 101 bp):

| Stage | Filename |
|-------|----------|
| Lane 1 merged | `{SAMPLE}_L1.fq.gz` |
| Lane 2 merged | `{SAMPLE}_L2.fq.gz` |
| Combined lanes | `{SAMPLE}.fq.gz`, `{SAMPLE}.fq` |
| FASTA | `{SAMPLE}.fa` |
| Length-normalized FASTA | `{SAMPLE}_{N}.fa` |
| 40-mers (plain text) | `{SAMPLE}_{N}.seq` |

Adapter FASTA is supplied separately (e.g. `Adapters_TruSeq3PE.fa`).

## Dependencies

- [`fastp`](https://github.com/OpenGene/fastp)
- [`seqtk`](https://github.com/lh3/seqtk)
- `gzip`, `bash`, `python3` (≥ 3.9 recommended)
- [Biopython](https://biopython.org/) (`pip install -r requirements.txt`)

## Usage

```bash
./preprocessing/run_preprocess.sh \
  --sample SAMPLE \
  --adapter-fasta ./Adapters_TruSeq3PE.fa \
  --ligand lig147 \
  --input-dir ./raw \
  --output-dir ./processed
```

Provide exactly one of `--ligand` or `--target-length` (no hard-coded default such as 154 or 101). Optionally widen the retention window and control trimming/padding:

```bash
./preprocessing/run_preprocess.sh \
  --sample SAMPLE \
  --adapter-fasta ./Adapters_TruSeq3PE.fa \
  --ligand lig200 \
  --min-input-len 145 \
  --max-input-len 163 \
  --trim-from right \
  --seed 1 \
  --input-dir ./raw \
  --output-dir ./processed
```

Equivalent standalone steps (mirroring the original workflow with generic names):

```bash
fastp -q 10 \
  --in1 SAMPLE_L1_R1.fastq.gz --in2 SAMPLE_L1_R2.fastq.gz \
  -m --merged_out SAMPLE_L1.fq.gz \
  --adapter_fasta=./Adapters_TruSeq3PE.fa \
  --json=./SAMPLE_L1.json --html=./SAMPLE_L1.html \
  --overlap_len_require=5 -gx --length_required=15 \
  --n_base_limit=10 -y --complexity_threshold=5

fastp -q 10 \
  --in1 SAMPLE_L2_R1.fastq.gz --in2 SAMPLE_L2_R2.fastq.gz \
  -m --merged_out SAMPLE_L2.fq.gz \
  --adapter_fasta=./Adapters_TruSeq3PE.fa \
  --json=./SAMPLE_L2.json --html=./SAMPLE_L2.html \
  --overlap_len_require=5 -gx --length_required=15 \
  --n_base_limit=10 -y --complexity_threshold=5

cat SAMPLE_L1.fq.gz SAMPLE_L2.fq.gz > SAMPLE.fq.gz
gzip -d -f SAMPLE.fq.gz
seqtk seq -a SAMPLE.fq > SAMPLE.fa
python3 preprocessing/filter_by_length.py SAMPLE.fa SAMPLE_101.fa --ligand lig147
python3 preprocessing/trim_to_kmer.py SAMPLE_101.fa SAMPLE_101.seq
```

### Helper scripts

`filter_by_length.py` is the main length-normalization helper. Provide either `--target_length` or `--ligand` (mutually exclusive; no hard-coded default). If `--min_input_len` / `--max_input_len` are omitted, both default to the resolved target length (exact-length retention only). Shorter reads in the window are padded with random A/T/C/G; longer reads are trimmed (`--trim_from right|left|center`).

```bash
# Normalize via ligand design key
python3 preprocessing/filter_by_length.py input.fa output_101.fa --ligand lig147

# Or set an explicit length
python3 preprocessing/filter_by_length.py input.fa output_154.fa --target_length 154

# Windowed retention + pad/trim (lig200 example)
python3 preprocessing/filter_by_length.py input.fa output_154.fa \
  --ligand lig200 \
  --min_input_len 145 \
  --max_input_len 163 \
  --trim_from right \
  --seed 1

# Extract 40-mers starting at offset 0 (default); use --center to center the window
python3 preprocessing/trim_to_kmer.py input_101.fa output.seq --kmer-size 40
```

## `fastp` parameters (current defaults)

| Flag | Value | Role |
|------|-------|------|
| `-q` | 10 | Base quality threshold |
| `-m` / `--merged_out` | — | Merge overlapping PE reads |
| `--adapter_fasta` | user-supplied | Adapter sequences |
| `--overlap_len_require` | 5 | Minimum overlap for merging |
| `-g` / `-x` | on | PolyG trimming / length filtering helpers as configured by fastp |
| `--length_required` | 15 | Discard very short merges |
| `--n_base_limit` | 10 | Max ambiguous bases |
| `-y` / `--complexity_threshold` | 5 | Low-complexity filter |

These defaults match the established SELEX preprocessing recipe used for this project and can be parameterised further as the pipeline evolves.
