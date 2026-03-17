#!/usr/bin/env python3
"""Detect 3' adapter sequence in small RNA-seq FASTQ files.

Searches the first N reads for common adapter sequences and reports which
adapter is present, along with the recommended --adapter_sequence setting.

Usage:
    detect_adapter.py sample.fastq.gz
    detect_adapter.py sample.fastq.gz --reads 20000
"""

import argparse
import gzip
import sys

ADAPTERS = {
    "Illumina TruSeq Universal": "AGATCGGAAGAGC",
    "TruSeq Small RNA (RA3)": "TGGAATTCTCGGGTGCCAAGG",
    "NEBNext Small RNA": "AGATCGGAAGAGCACACGTCT",
    "Nextera": "CTGTCTCTTATACACATCT",
}

# Corresponding --adapter_sequence recommendations
RECOMMENDATIONS = {
    "Illumina TruSeq Universal": None,  # fastp auto-detect handles this
    "TruSeq Small RNA (RA3)": "TGGAATTCTCGGGTGCCAAGG",
    "NEBNext Small RNA": None,  # fastp auto-detect handles this
    "Nextera": "CTGTCTCTTATACACATCT",
}


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def detect(fastq_path, n_reads=10000):
    counts = {name: 0 for name in ADAPTERS}
    total = 0

    with open_fastq(fastq_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:  # sequence lines only
                continue
            seq = line.strip()
            total += 1
            for name, motif in ADAPTERS.items():
                if motif in seq:
                    counts[name] += 1
            if total >= n_reads:
                break

    return counts, total


def main():
    parser = argparse.ArgumentParser(
        description="Detect 3' adapter in small RNA-seq FASTQ"
    )
    parser.add_argument("fastq", help="FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument(
        "--reads", type=int, default=10000, help="Number of reads to scan (default: 10000)"
    )
    args = parser.parse_args()

    print(f"Scanning first {args.reads} reads of {args.fastq}...\n")
    counts, total = detect(args.fastq, args.reads)

    if total == 0:
        print("ERROR: No reads found in file.", file=sys.stderr)
        sys.exit(1)

    print(f"Reads scanned: {total:,}\n")
    print(f"{'Adapter':<30} {'Matches':>8} {'Rate':>8}")
    print("-" * 50)
    for name, count in sorted(counts.items(), key=lambda x: -x[1]):
        rate = 100 * count / total
        print(f"{name:<30} {count:>8,} {rate:>7.1f}%")

    # Find best match
    best_name = max(counts, key=counts.get)
    best_count = counts[best_name]
    best_rate = 100 * best_count / total

    print()
    if best_rate < 5:
        # Check if reads have variable lengths (pre-trimmed)
        print("No adapter detected (< 5% match rate).")
        print("Data may be pre-trimmed — use the default (--adapter_sequence not needed).")
    else:
        rec = RECOMMENDATIONS[best_name]
        print(f"Detected: {best_name} ({best_rate:.0f}% of reads)")
        if rec is None:
            print("Recommendation: use the default (fastp auto-detect handles this adapter)")
        else:
            print(f"Recommendation: --adapter_sequence '{rec}'")
            print("  (fastp cannot auto-detect this adapter for single-end reads)")


if __name__ == "__main__":
    main()
