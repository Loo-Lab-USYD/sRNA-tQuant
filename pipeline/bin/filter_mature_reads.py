#!/usr/bin/env python3
"""Replace removePrecursor.pl: keep only reads mapping to mature tRNA region.

Usage: filter_mature_reads.py <pre_trnas.bed12> <input.sam> <input.fastq> [--threads N] [--flank F]

Writes qualifying FASTQ records to stdout (matching legacy behaviour).

Logic (matching legacy removePrecursor.pl exactly):
  - Load pre-tRNA BED12: build key as 'name(strand)', value is sum of block sizes
  - For each SAM alignment:
    - Parse CIGAR to compute insertions (I) and deletions (D)
    - end = start + read_length - 1 + D - I
    - If start > flank AND end <= tRNA_length - flankEnd: read maps to mature region
    - flankEnd = flank - 6 (accounts for CCACCA tail)
  - Extract qualifying reads from FASTQ and write to stdout
"""

import argparse
import gzip
import os
import re
from concurrent.futures import ProcessPoolExecutor


def parse_cigar_indels(cigar_str):
    """Parse CIGAR string, return total insertions and deletions."""
    parts = re.split(r"([IMXD])", cigar_str)
    insertions = 0
    deletions = 0
    for i, part in enumerate(parts):
        if part == "I":
            insertions += int(parts[i - 1])
        elif part == "D":
            deletions += int(parts[i - 1])
    return insertions, deletions


def _scan_sam_chunk(chunk_args):
    """Scan a byte range of the SAM file for mature-region reads."""
    sam_path, start_offset, end_offset, trna_length, flank, flank_end = chunk_args
    mature = set()
    with open(sam_path, "rb") as fh:
        fh.seek(start_offset)
        if start_offset > 0:
            fh.readline()

        while fh.tell() < end_offset:
            line = fh.readline()
            if not line:
                break
            if line.startswith(b"@"):
                continue
            fields = line.strip().split(b"\t")
            if len(fields) < 10:
                continue
            read_name = fields[0].decode("ascii")
            ref_name = fields[2].decode("ascii")
            start = int(fields[3])
            seq = fields[9]
            cigar = fields[5].decode("ascii")

            ins, dels = parse_cigar_indels(cigar)
            end = start + len(seq) - 1 + dels - ins

            # SAM ref_name may be full FASTA header:
            #   "chr1.tRNA149-UndetNNN::chr1:7930228-7930398(-)"
            # BED12 key is "chr1.tRNA149-UndetNNN(-)"
            if "::" in ref_name:
                short_name = ref_name.split("::")[0]
                strand = ref_name.rstrip(")").rsplit("(", 1)[-1]
                lookup_key = f"{short_name}({strand})"
            else:
                lookup_key = ref_name

            if lookup_key in trna_length:
                if start > flank and end <= trna_length[lookup_key] - flank_end:
                    mature.add(read_name)

    return mature


def filter_mature_reads(bed12_path, sam_path, fastq_path, flank=50, threads=1):
    flank_end = flank - 6  # CCACCA tail

    # Load pre-tRNA lengths from BED12
    trna_length = {}
    opener = gzip.open if bed12_path.endswith(".gz") else open
    with opener(bed12_path, "rt") as fh:
        for line in fh:
            fields = line.strip().split("\t")
            name = fields[3]
            strand = fields[5]
            key = f"{name}({strand})"
            block_sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
            trna_length[key] = sum(block_sizes)

    if threads <= 1:
        # Single-threaded: original logic
        mature_reads = set()
        opener = gzip.open if sam_path.endswith(".gz") else open
        with opener(sam_path, "rt") as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                fields = line.strip().split("\t")
                read_name = fields[0]
                ref_name = fields[2]
                start = int(fields[3])
                seq = fields[9]
                cigar = fields[5]

                ins, dels = parse_cigar_indels(cigar)
                end = start + len(seq) - 1 + dels - ins

                # SAM ref_name may be full FASTA header
                if "::" in ref_name:
                    short_name = ref_name.split("::")[0]
                    strand = ref_name.rstrip(")").rsplit("(", 1)[-1]
                    lookup_key = f"{short_name}({strand})"
                else:
                    lookup_key = ref_name

                if lookup_key in trna_length:
                    if start > flank and end <= trna_length[lookup_key] - flank_end:
                        mature_reads.add(read_name)
    else:
        # Parallel SAM scanning
        file_size = os.path.getsize(sam_path)
        chunk_size = max(1, file_size // threads)
        chunks = []
        for i in range(threads):
            start_off = i * chunk_size
            end_off = file_size if i == threads - 1 else (i + 1) * chunk_size
            chunks.append((sam_path, start_off, end_off, trna_length, flank, flank_end))

        mature_reads = set()
        with ProcessPoolExecutor(max_workers=threads) as executor:
            for chunk_mature in executor.map(_scan_sam_chunk, chunks):
                mature_reads.update(chunk_mature)

    # Extract qualifying reads from FASTQ
    fq_opener = gzip.open if fastq_path.endswith(".gz") else open
    with fq_opener(fastq_path, "rt") as fh:
        entry = []
        for line in fh:
            line = line.rstrip("\n")
            entry.append(line)
            if len(entry) == 4:
                header, seq, plus_line, qual = entry
                entry = []
                # Extract read name: strip '@' prefix and description after first space
                read_id = header[1:].split()[0]
                if read_id in mature_reads:
                    print(f"{header}\n{seq}\n{plus_line}\n{qual}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter reads to mature tRNA region")
    parser.add_argument("bed12", help="Pre-tRNA BED12 file")
    parser.add_argument("sam", help="Input SAM file")
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument("--flank", type=int, default=50, help="Flank length (default: 50)")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel workers")
    args = parser.parse_args()

    filter_mature_reads(args.bed12, args.sam, args.fastq, args.flank, args.threads)
