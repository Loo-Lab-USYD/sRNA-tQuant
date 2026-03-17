# Applying the tRNA Pipeline to New Datasets

Practical guide for running the pipeline on public or in-house small RNA-seq data. Covers adapter detection, parameter choices, and known dataset-specific quirks.

## Quick Start Checklist

Before running any new dataset:

1. **Identify the library prep** (GEO "Library strategy" / "Construction protocol")
2. **Check the adapter** (see [Adapter Detection](#adapter-detection) below)
3. **Check the read length** (determines `--max_length_se` / `--max_length_pe`)
4. **Run on a small subset first** (subsample 50K reads, use mini genome)

## Compatible Data Types

### What works

The pipeline quantifies **mature cytoplasmic tRNAs at anticodon resolution**. It works with any sequencing data that captures tRNA-derived fragments in the ~18-50nt range:

| Data Type | Compatible? | Notes |
|---|---|---|
| **Small RNA-seq** (TruSeq, NEBNext) | Yes | Primary target. Captures tRNA fragments 18-35nt. |
| **Hydro-tRNAseq** | Yes | Alkaline hydrolysis of full-length tRNAs → 35-50nt fragments. May need `--max_length_se 60`. |
| **TGIRT-seq / mim-tRNAseq** | No | Uses specialised RT that reads through modifications. Different quantification model — use their own pipeline. |
| **YAMAT-seq** | Possibly | Y-shaped adapter captures full-length tRNAs. Untested — would need longer `--max_length_se`. |
| **Total RNA-seq** | Poor | tRNA fragments are a tiny fraction. Low sensitivity, but technically works if you have deep sequencing. |
| **CLIP-seq / Ribo-seq** | No | Different biology (protein-bound RNA / ribosome footprints). |

### Insert size constraints

The pipeline applies a post-trimming length filter that determines which reads are quantified:

```
--min_length 10      # Discard reads < 10bp after trimming (too short to map)
--max_length_se 50   # Discard SE reads > 50bp after trimming
--max_length_pe 125  # Discard PE reads > 125bp after trimming
```

These defaults are tuned for **standard small RNA-seq** where tRNA-derived fragments are 18-35nt. Here's why the bounds matter:

**Too short (< 10bp)**: Can't be uniquely mapped. These are adapter dimers or degradation artifacts.

**Too long (> 50bp for SE)**: Mature tRNAs are 73-93nt. Standard small RNA-seq captures fragments, not full-length tRNAs — most inserts are 18-35nt from RT truncation at modification sites. Reads that are still > 50bp after adapter trimming are likely:
- Non-tRNA reads (rRNA, mRNA fragments, miRNA precursors)
- Full-length tRNAs (rare in standard small RNA-seq, common in Hydro-tRNAseq/YAMAT-seq)

**When to adjust `--max_length_se`**:
- **Hydro-tRNAseq**: Produces 35-50nt fragments from alkaline hydrolysis. Set `--max_length_se 60`.
- **YAMAT-seq / full-length capture**: Set `--max_length_se 100` to retain full-length tRNA reads (~75-93nt).
- **Don't increase for standard small RNA-seq**: The 50bp filter correctly removes non-tRNA contaminants. Making it larger adds noise without signal.

### What about the anticodon-level vs transcript-level distinction?

The pipeline clusters tRNA genes by anticodon (e.g., all tRNA-Ala-AGC genes → one cluster). This is appropriate because:
- Standard small RNA-seq reads are truncated at modification sites, producing fragments too short to distinguish individual tRNA genes within an anticodon family
- Anticodon-level resolution is what matters for codon usage studies

The bowtie2 seed parameters control resolution:
```
--bowtie2_seed_length 14        # Anticodon-level (default, less stringent)
--bowtie2_seed_length 18        # Transcript-level (more stringent, loses some reads)
--bowtie2_seed_extensions 25    # Default; increase to 50 for transcript-level
```

Transcript-level quantification requires longer reads and/or TGIRT-based library prep that reads through modifications.

## Adapter Detection

This is the single most important step. Getting the adapter wrong means fastp can't trim reads, untrimmed reads won't align, and you'll quantify essentially nothing with no error message.

### The Problem

Small RNA-seq library preps ligate RNA adapters to the small RNA insert, then add Illumina sequencing adapters during PCR. The adapter that appears in the **sequencing reads** depends on the library prep:

| Library Prep | Adapter in Reads | fastp Auto-Detect? | `--adapter_sequence` |
|---|---|---|---|
| Most standard small RNA-seq | Illumina universal: `AGATCGGAAGAGC...` | Yes (SE and PE) | `null` (default) |
| TruSeq Small RNA (some datasets) | RA3: `TGGAATTCTCGGGTGCCAAGG` | PE only, **not SE** | Must pass explicitly |
| NEBNext Small RNA | `AGATCGGAAGAGCACACGTCT` | Yes | `null` |
| Unknown | Check reads first | -- | -- |

### How to Check

Inspect the first few reads for adapter sequence:

```bash
# Look at raw reads — adapter appears after the insert
zcat sample.fastq.gz | head -40 | awk 'NR%4==2'
```

Common patterns:
```
# Illumina universal adapter (fastp auto-detect works):
TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTAGATCGGAAGAGCACAC
                                  ^^^^^^^^^^^^^^^^
# TruSeq Small RNA RA3 adapter (must pass explicitly):
TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTTGGAATTCTCGGGTGCCA
                                  ^^^^^^^^^^^^^^^^^^^^^
```

You can also count adapter occurrences:

```bash
# Check for Illumina universal
zcat sample.fastq.gz | head -2000 | grep -c "AGATCGGAAGAG"

# Check for TruSeq Small RNA RA3
zcat sample.fastq.gz | head -2000 | grep -c "TGGAATTCTCGGGTGCCAAGG"
```

If neither is found, the data may be pre-trimmed (check if reads vary in length).

### Pipeline Usage

```bash
# Auto-detect (works for Illumina universal adapter) — the default
nextflow run pipeline/main.nf --input samplesheet.csv ...

# Explicit adapter (required for TruSeq Small RNA RA3 on SE reads)
nextflow run pipeline/main.nf --input samplesheet.csv \
    --adapter_sequence 'TGGAATTCTCGGGTGCCAAGG' ...
```

### Why fastp Can't Auto-Detect RA3 for SE Reads

fastp's SE auto-detection algorithm looks for the Illumina universal adapter (`AGATCGGAAGAG`) by analyzing overlap patterns. The TruSeq Small RNA RA3 adapter (`TGGAATTCTCGGGTGCCAAGG`) is a different sequence that the algorithm doesn't search for. For PE reads, fastp uses insert-size analysis which is adapter-agnostic, so it works regardless.

## Tested Datasets

### GSE137834 — Hernandez-Alias et al. (5 cell lines, small RNA-seq + Hydro-tRNAseq)

- **Platform**: HiSeq 2500, 50bp SE
- **Adapter**: Illumina universal (`AGATCGGAAGAG`) — fastp auto-detect works
- **`--adapter_sequence`**: `null` (default)
- **Samples**: HEK293, HeLa, HCT116, BJ, MDA-MB-231 (3 reps each)
- **Notes**:
  - ~98% of reads contain adapter (short inserts in 50bp reads)
  - ~5% of reads are tRNA-derived (~900K per sample out of ~18M)
  - Both smallRNAseq and hydroseq columns in the published matrix

### ENCODE ENCSR000AES — K562 small RNA-seq

- **Platform**: HiSeq 2000, 101bp SE
- **Adapter**: TruSeq Small RNA RA3 (`TGGAATTCTCGGGTGCCAAGG`) — must pass explicitly
- **`--adapter_sequence`**: `TGGAATTCTCGGGTGCCAAGG`
- **Samples**: K562 (2 reps: ENCFF001REQ, ENCFF001REL)
- **Notes**:
  - TAP-treated, rRNA-depleted total RNA <200nt
  - 101bp reads mean more adapter readthrough — ~78% of reads contain adapter
  - After trimming, reads are typically 18-35bp (standard small RNA insert sizes)
  - `--max_length_se 50` in the default config filters appropriately after trimming

### GSE152621 / mim-tRNAseq — Behrens et al. (HEK293T, K562, iPSC)

- **Used as**: Gold standard reference for benchmarking (not run through our pipeline)
- **Method**: TGIRT reverse transcriptase (reads through modifications)
- **Anticodon count table**: `GSE152621_Hsap_anticodon_counts.csv`
- **Naming**: `Ala-AGC` format — harmonise by removing hyphen to match our `AlaAGC`
- **Notes**:
  - Includes mitochondrial tRNAs (prefix `mito-`) and eColiLys spike-in — filter these
  - Raw counts, not RPM — normalise before comparison
  - Cross-study correlation with standard small RNA-seq is typically rho ~ 0.5-0.7 due to lab/normalization effects, not method failure

## Read Length Considerations

| Read Length | Typical Source | Adapter Situation | Notes |
|---|---|---|---|
| 36bp | Older HiSeq | Most reads have adapter | Very short inserts only |
| 50bp | HiSeq 2500 | ~95%+ reads have adapter | Most common for small RNA-seq |
| 75bp | NextSeq/NovaSeq | ~80%+ reads have adapter | Good for capturing full tRNA fragments |
| 101bp | HiSeq 2000/ENCODE | ~75%+ reads have adapter | Lots of adapter readthrough |
| 150bp | NovaSeq | ~70%+ reads have adapter | Excessive for small RNA |

The pipeline's default `--max_length_se 50` filters reads after trimming. If your inserts are genuinely longer (e.g., Hydro-tRNAseq produces ~35-50bp fragments), you may need to increase this.

## Pre-Built vs De Novo Reference

```bash
# First run: build reference from genome FASTA (~40 min for hg38)
nextflow run pipeline/main.nf --genome_fasta hg38.fa --outdir results ...

# The reference bundle is automatically published to results/reference/
# Subsequent runs: reuse it to skip the 40-minute build
nextflow run pipeline/main.nf --trna_reference results/reference/ ...
```

The reference bundle is genome-specific. An hg38 bundle works for any hg38 dataset regardless of library prep or read length.

## Common Pitfalls

### 1. Silent adapter failure

**Symptom**: Pipeline completes successfully but READS count in RPM matrix is very low (< 1000 for a typical dataset).

**Cause**: Wrong adapter or no adapter trimming. Reads retain adapter sequence, can't align to tRNA reference.

**Check**: Look at the fastp JSON in the work directory:
```bash
python3 -c "
import json
with open('work/<hash>/<sample>_fastp.json') as f:
    d = json.load(f)
print('Adapter-trimmed reads:', d['adapter_cutting']['adapter_trimmed_reads'])
print('Total reads:', d['summary']['before_filtering']['total_reads'])
"
```
If adapter-trimmed reads is < 1% of total, the adapter is wrong.

### 2. Pre-trimmed data from SRA/GEO

Some submitters trim adapters before uploading to SRA. Signs:
- Reads have variable lengths (not all the same)
- No adapter sequence found in reads
- fastp adapter_trimmed_reads is ~0

This is fine — the pipeline handles it. Just use `--adapter_sequence null` (default).

### 3. Modification calling sample mismatch

Fixed in the current version. The pipeline now joins channels by sample ID before passing to ADJUST_MODIFICATIONS and COUNT_ALLELES. If you see a sample's RPM being adjusted with another sample's ASE data in older versions, update the pipeline.

### 4. Nextflow cache after parameter changes

Changing `--adapter_sequence` changes the TRIM_READS task hash, so `-resume` correctly re-runs trimming. But if you change the adapter in `nextflow.config` (not on the command line), you need to verify the cache was actually busted:
```bash
# Check the TRIM_READS command that actually ran:
cat work/<hash>/.command.sh | grep adapter
```

## Samplesheet Format

```csv
sample_id,fastq_1,fastq_2
HEK293_rep1,/absolute/path/to/sample_R1.fastq.gz,
HEK293_rep2,/absolute/path/to/sample_R1.fastq.gz,
```

- Paths must be absolute
- Leave `fastq_2` empty for SE data (keep the trailing comma)
- `sample_id` becomes the column name in the output RPM matrix
