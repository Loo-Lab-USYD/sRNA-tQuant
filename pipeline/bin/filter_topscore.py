#!/usr/bin/env python3
"""Filter SAM/BAM to retain only equal top-scoring alignments per read.

Usage: filter_topscore.py <input.sam> <output.sam>

Required for the Salmon quantification path. When Bowtie2 is run with -a
(report all valid alignments), it outputs suboptimal alignments alongside
the best ones. Salmon expects only the best-scoring set per read.

Logic:
  - Group alignments by read name (input must be name-sorted or name-grouped)
  - For each read, find the maximum alignment score (AS:i tag)
  - Output only alignments with score equal to the maximum
  - SAM headers are passed through unchanged

Reference: Smith et al. (2024) — "The alignments must then be filtered to
retain only the equal top-scoring alignments for each read, before passing
the alignments to salmon quant."
"""

import re
import sys


def get_alignment_score(sam_line):
    """Extract AS:i:<score> from SAM optional fields."""
    match = re.search(r"AS:i:(-?\d+)", sam_line)
    if match:
        return int(match.group(1))
    return None


def filter_topscore(sam_in, sam_out):
    with open(sam_in) as in_fh, open(sam_out, "w") as out_fh:
        current_read = None
        alignments = []
        scores = []

        def flush_group():
            """Write only alignments with the best score."""
            if not alignments:
                return
            max_score = max(scores)
            for aln, score in zip(alignments, scores):
                if score == max_score:
                    out_fh.write(aln)

        for line in in_fh:
            if line[0] == "@":
                out_fh.write(line)
                continue

            fields = line.split("\t")
            read_name = fields[0]
            score = get_alignment_score(line)

            if read_name == current_read:
                alignments.append(line)
                scores.append(score if score is not None else 0)
            else:
                flush_group()
                current_read = read_name
                alignments = [line]
                scores = [score if score is not None else 0]

        # Flush final group
        flush_group()


if __name__ == "__main__":
    filter_topscore(sys.argv[1], sys.argv[2])
