# CLAUDE.md

The role of this file is to describe common mistakes and confusion points that agents might encounter as they work in this project. If you ever encounter something in the project that surprises you, please alert the developer working with you and indicate that his is the case in the CLAUDE.md file to help prevent future agents from having the same issue.

## Known Gotchas

### rRNA depletion (REMOVE_RRNA)
The optional `--rrna_index` step removes rRNA reads with `bowtie2 --very-fast --un-gz` before mapping to the artificial genome. This is purely a noise-reduction step — rRNA reads can't map to tRNA anyway, but removing them upfront can reduce Bowtie2 runtime. Provide a bowtie2 index built from rRNA reference sequences. Use the `rrna_index` directory (must contain `rrna.1.bt2`, etc.). When absent (default), the step is skipped entirely.

### Fastp max_length_se
`max_length_se` was increased from 50 → 75 (2026-05-02) because many modern sRNA-seq datasets produce longer reads (e.g. 75bp PE). The old 50bp cutoff was too conservative and could discard valid full-length tRNA reads. For very old datasets (<50bp reads), this has no effect since fastp already trims to read length.

### fastp adapter detection — DO NOT assume the adapter type
fastp's auto-adapter detection works well for the Illumina TruSeq universal adapter (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`), which is the most common adapter in small RNA-seq datasets (including GSE137834). Passing the wrong adapter (e.g. TruSeq Small RNA `TGGAATTCTCGGGTGCCAAGG` on data that uses the universal adapter) causes fastp to skip trimming, leaving reads at full length with adapter contamination. This makes them unmappable to 73-77 bp tRNA clusters, resulting in <1% alignment rate and garbage RPMs.

**The safe default is `adapter_sequence = null`** (fastp auto-detect). Only override with `--adapter_sequence` if you have confirmed the adapter from the library prep protocol or from the fastp auto-detect report. The pipeline config default is already `null`.

### Salmon minAssignedFrags
Salmon's default `--minAssignedFrags 10` will fail on small test datasets. The pipeline config now sets `--minAssignedFrags 1` to handle edge cases gracefully.

### adjust_modifications.py zero-depth positions
The ASE (allele-specific expression) data can have entries with `rawDepth == 0`. The `adjust_modifications.py` script must guard against division by zero when computing `altCount / rawDepth`.

### Unit test Python environment
Unit tests that call `subprocess.run(['python3', ...])` will fail if system Python lacks pandas. Use the `trna_validation` mamba env: `/home/labhund/mamba_envs/trna_validation/bin/python`.

### Nextflow cache invalidation
Bash comments don't change Nextflow's task hash. To bust cache after changing script logic, use `echo "vN" > /dev/null` with an incrementing version number inside the process script block.
