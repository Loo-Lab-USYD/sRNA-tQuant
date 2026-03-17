#!/usr/bin/env python3
"""Replace addCCA.pl: append CCA tail to mature tRNAs, filter pseudogenes.

Usage: add_cca.py <input.fa> <output.fa>

Logic (matching legacy addCCA.pl exactly):
  - Read multi-line FASTA
  - Skip entries whose header contains 'pseudo' (case-insensitive)
  - Lowercase the sequence and append 'cca'
  - Write single-line FASTA preserving original entry order
"""

import sys


def add_cca(in_path, out_path):
    headers = []
    seqs = {}

    with open(in_path) as fh:
        current_id = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current_id = line
                headers.append(current_id)
                seqs[current_id] = ""
            else:
                seqs[current_id] += line

    with open(out_path, "w") as fh:
        for header in headers:
            if "pseudo" not in header.lower():
                fh.write(f"{header}\n")
                fh.write(f"{seqs[header].lower()}cca\n")


if __name__ == "__main__":
    add_cca(sys.argv[1], sys.argv[2])
