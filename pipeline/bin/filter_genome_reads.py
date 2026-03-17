#!/usr/bin/env python3
"""Replace removeGenomeMapper.pl: discard reads that map to the genome.

Usage: filter_genome_reads.py <pre_trnas.fa> <input.sam> <output.sam> [--threads N] [--flank F]

Logic (matching legacy removeGenomeMapper.pl exactly):
  - Load pre-tRNA contig names from FASTA headers
  - For each SAM alignment:
    - If the read sequence length >= 30 (or >= flank):
      - If the alignment target is NOT a pre-tRNA contig, flag the read
        for removal (it mapped to the genome)
    - Reads shorter than 30 are never flagged
  - A read is discarded if ANY of its alignments hit a non-pre-tRNA contig
  - Output all SAM lines for reads that were never flagged

Streaming two-pass: first pass only collects flagged read names (not full
SAM lines) to avoid loading multi-GB SAM files into memory.
"""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor


def _scan_chunk(chunk_args):
    """Scan a byte range of the SAM file for genome-mapping reads."""
    sam_path, start_offset, end_offset, pretrna_ids, min_len = chunk_args
    discard = set()
    with open(sam_path, "rb") as fh:
        # Seek to start; if not at beginning, skip partial line
        fh.seek(start_offset)
        if start_offset > 0:
            fh.readline()

        while fh.tell() < end_offset:
            line = fh.readline()
            if not line:
                break
            if line.startswith(b"@"):
                continue
            fields = line.split(b"\t", 11)
            if len(fields) < 10:
                continue
            read_name = fields[0].decode("ascii")
            ref_name = fields[2].decode("ascii")
            seq = fields[9]

            if len(seq) >= min_len and ref_name not in pretrna_ids:
                discard.add(read_name)

    return discard


def filter_genome_reads(pretrna_fa, sam_in, sam_out, flank=50, threads=1):
    # Load pre-tRNA contig names
    pretrna_ids = set()
    with open(pretrna_fa) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                pretrna_ids.add(line[1:])

    min_len = min(flank, 30)  # effectively 30 since flank=50

    if threads <= 1:
        # Single-threaded: original logic
        discard_reads = set()
        with open(sam_in) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                fields = line.split("\t", 11)  # only split enough fields
                read_name = fields[0]
                ref_name = fields[2]
                seq = fields[9]

                if len(seq) >= min_len and ref_name not in pretrna_ids:
                    discard_reads.add(read_name)
    else:
        # Parallel first pass: chunk file by byte offset
        file_size = os.path.getsize(sam_in)
        chunk_size = max(1, file_size // threads)
        chunks = []
        for i in range(threads):
            start = i * chunk_size
            end = file_size if i == threads - 1 else (i + 1) * chunk_size
            chunks.append((sam_in, start, end, pretrna_ids, min_len))

        discard_reads = set()
        with ProcessPoolExecutor(max_workers=threads) as executor:
            for chunk_discard in executor.map(_scan_chunk, chunks):
                discard_reads.update(chunk_discard)

    # Second pass: stream through SAM again, output kept reads
    with open(sam_in) as fh_in, open(sam_out, "w") as fh_out:
        for line in fh_in:
            if line.startswith("@"):
                fh_out.write(line)
                continue
            read_name = line.split("\t", 1)[0]
            if read_name not in discard_reads:
                fh_out.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter genome-mapping reads from SAM")
    parser.add_argument("pretrna_fa", help="Pre-tRNA FASTA file")
    parser.add_argument("sam_in", help="Input SAM file")
    parser.add_argument("sam_out", help="Output SAM file")
    parser.add_argument("--flank", type=int, default=50, help="Flank length (default: 50)")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel workers")
    args = parser.parse_args()

    filter_genome_reads(args.pretrna_fa, args.sam_in, args.sam_out, args.flank, args.threads)
