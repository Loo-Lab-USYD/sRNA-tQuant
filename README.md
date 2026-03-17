# tRNA Quantification from Small RNA-seq

Quantify tRNA abundances at anticodon resolution from standard small RNA-seq data. Designed for generating cell-line-specific tRNA profiles for codon optimisation and translational research.

## Quick start

```bash
# Build container
cd pipeline && docker build -t trna-mapping:latest . && cd ..

# Run
nextflow run pipeline/main.nf \
  --input samplesheet.csv \
  --genome_fasta /path/to/hg38.fa \
  --outdir results \
  -profile docker -resume
```

Output: `results/results_rpm.csv` — anticodon × sample RPM matrix.

## Samplesheet

```csv
sample_id,fastq_1,fastq_2
HeLa_rep1,/data/HeLa_rep1.fastq.gz,
K562_rep1,/data/K562_rep1_R1.fastq.gz,/data/K562_rep1_R2.fastq.gz
```

`fastq_2` is optional (leave empty for single-end).

## Adapter detection

Getting the adapter right is critical. The wrong adapter means fastp can't trim reads, untrimmed reads won't align, and the pipeline produces near-zero RPMs with no error message. The pipeline will warn you if < 2% of reads are identified as tRNA, but it's better to check upfront.

**Check your adapter before running:**

```bash
python3 pipeline/bin/detect_adapter.py sample.fastq.gz
```

**Common library preps and their settings:**

| Library Prep | Adapter | `--adapter_sequence` | fastp Auto-detect? |
|---|---|---|---|
| TruSeq (most standard small RNA-seq) | Illumina universal | `null` (default) | Yes |
| NEBNext Small RNA | Illumina-compatible | `null` (default) | Yes |
| TruSeq Small RNA (some ENCODE data) | RA3: `TGGAATTCTCGGGTGCCAAGG` | `'TGGAATTCTCGGGTGCCAAGG'` | PE only, **not SE** |
| Pre-trimmed (SRA/GEO) | None in reads | `null` (default) | N/A |

For most datasets, the default (`null` = fastp auto-detect) works. Only override if you have confirmed the adapter or `detect_adapter.py` tells you to.

## Getting started

### 1. Install prerequisites

The pipeline requires [Nextflow](https://www.nextflow.io/docs/latest/install.html) (>= 23.04) and Java 11+. All bioinformatics tools (bowtie2, Salmon, tRNAscan-SE, etc.) run inside a container or conda environment — you don't need to install them separately.

Choose one execution environment:

**Docker** (recommended):
```bash
cd pipeline && docker build -t trna-mapping:latest . && cd ..
```

**Conda/Mamba** (no Docker needed):
```bash
# The pipeline creates the environment automatically from pipeline/environment.yml
# Just use -profile conda when running
```

**Singularity** (HPC clusters):
```bash
singularity build trna-mapping.sif docker-daemon://trna-mapping:latest
```

### 2. Download a reference genome

The pipeline needs a reference genome FASTA to predict and extract tRNA sequences. For human (hg38):

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

Other genomes are available from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) or [Ensembl](https://www.ensembl.org/info/data/ftp/index.html). Any soft-masked or unmasked FASTA works — the pipeline runs its own masking.

### 3. First run — build the reference (~40 min)

The first run builds a tRNA reference bundle automatically from the genome FASTA using tRNAscan-SE, then indexes it with bowtie2. This takes ~40 minutes for hg38 and is cached by Nextflow's `-resume`.

```bash
nextflow run pipeline/main.nf \
  --input samplesheet.csv \
  --genome_fasta hg38.fa \
  --outdir results \
  -profile docker -resume
```

### 4. Reuse the reference bundle for future runs

The first run publishes the built reference to `results/reference/`. Subsequent runs can skip the 40-minute build by pointing to it:

```bash
nextflow run pipeline/main.nf \
  --input new_samplesheet.csv \
  --trna_reference results/reference/ \
  --outdir results \
  -profile docker -resume
```

The bundle is genome-specific — an hg38 bundle works for any hg38 dataset regardless of library prep or read length.

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quantifier` | `salmon` | `salmon` (EM-based, recommended) or `decision` (legacy) |
| `--span_weighted` | `false` | Blend Salmon with anticodon-spanning reads for improved isoacceptor resolution |
| `--skip_modification_calling` | `true` | Skip variant-based modification detection |
| `--adapter_sequence` | `null` | 3' adapter (null = fastp auto-detect, recommended) |
| `--max_cpus` | 8 | Maximum CPUs per process |
| `--max_memory` | `16.GB` | Maximum memory per process |

## Span-weighted quantification

Standard small RNA-seq captures tRNA fragments (median ~16 bp), most of which cannot distinguish isoacceptors of the same amino acid. The `--span_weighted` mode blends Salmon's EM estimates with direct counts from the ~7–10% of reads that span the anticodon position, improving isoacceptor resolution where it matters most.

Validated on HEK293, K562, and iPSC against mim-tRNAseq (TGIRT) gold standard:

| Cell line | Salmon ρ | Span-weighted ρ | Improvement |
|-----------|---------|-----------------|-------------|
| HEK293 | 0.697 | 0.746 | +0.049 |
| K562 | 0.542 | 0.594 | +0.053 |
| iPSC | 0.535 | 0.556 | +0.021 |

Validation methodology described in [Applying to New Datasets](docs/applying-to-new-datasets.md).

### Validation data

| Cell line | Source | Accessions | Reference (gold standard) |
|-----------|--------|-----------|--------------------------|
| HEK293 (3 reps) | [GSE137834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137834) | SRR10162601–603 | mim-tRNAseq HEK293T ([GSE152621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152621)) |
| K562 (2 reps) | [ENCODE](https://www.encodeproject.org/) | ENCFF001REQ, ENCFF001REL | mim-tRNAseq K562 (GSE152621) |
| iPSC (2 reps) | ENCODE | ENCFF004BLT, ENCFF450POH | mim-tRNAseq iPSC (GSE152621) |
| hESC (2 reps) | ENCODE | ENCFF990ZUE, ENCFF307UVF | Pending |

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Docker container |
| `singularity` | Singularity `.sif` |
| `gadi` | NCI Gadi cluster (PBS Pro + Singularity) |
| `conda` | Conda environment |
| `test` | Minimal test (2 CPUs, 4 GB) |

## Output

```
results/
├── results_rpm.csv            # Anticodon × sample RPM matrix
├── results_rpm_nomod.csv      # RPM without modification adjustment
├── modifications.csv          # Detected modifications (if enabled)
├── multiqc_report.html        # QC summary
├── reference/                 # Reusable tRNA reference bundle (from --genome_fasta runs)
└── pipeline_info/             # Timeline, resource report, DAG
```

## Pipeline overview

```
genome.fa ─► PREPARE_REFERENCE (tRNAscan-SE → mask → cluster → index)
                      │
sample.fq ─► TRIM_READS (fastp)
                      │
              MAP_ARTIFICIAL (bowtie2 → filter genome → filter mature)
                      │
              MAP_CLUSTERS (bowtie2 to clustered tRNA reference)
                      │
              ┌── Salmon path (default) ──────────────────┐
              │  Filter top-scoring → Salmon EM → RPM     │
              │  Optional: span-weighted blending          │
              ├── Decision path (legacy) ─────────────────┤
              │  Resolve multimappers by anticodon → RPM  │
              └───────────────────────────────────────────┘
                      │
              AGGREGATE_RESULTS → MultiQC
```

## Tests

Unit tests require Python 3 with pandas.

```bash
python3 -m pytest tests/ -v

# Integration tests (requires Nextflow + Docker)
bash tests/integration/run_integration_test.sh
bash tests/integration/run_k562_integration_test.sh
```

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) >= 23.04
- Java 11+
- One of: Docker, Singularity, or Conda/Mamba (see [Getting started](#getting-started))

## Acknowledgements

The core tRNA quantification strategy — mapping to an artificial genome, filtering by mature tRNA region, and clustering by anticodon — originates from [Anne Hoffmann's tRNA mapping pipeline](https://github.com/AnneHoffmann/tRNA_mapping) and was subsequently used by [Hernandez-Alias et al. (2020)](https://doi.org/10.15252/msb.20199896). This pipeline is a ground-up reimplementation in Nextflow/Python with probabilistic quantification via Salmon, span-weighted isoacceptor resolution, and automated reference building.

## References

- Hernandez-Alias X et al. "Translational efficiency across healthy and tumor tissues is determined by codon composition of tissue-specific human transcriptomes." *Mol Syst Biol.* 2020. https://doi.org/10.15252/msb.20199896
- Patro R et al. "Salmon provides fast and bias-aware quantification of transcript expression." *Nat Methods.* 2017. https://doi.org/10.1038/nmeth.4197
- Behrens A et al. "High-resolution quantitative profiling of tRNA abundance and modification status in eukaryotes by mim-tRNAseq." *Molecular Cell.* 2021. https://doi.org/10.1016/j.molcel.2021.01.028
- Smith KP et al. "Benchmarking tRNA-Seq quantification approaches exploiting a Saccharomyces cerevisiae dataset." *eLife.* 2024.

## License

[MIT](LICENSE) - Markus Williams
