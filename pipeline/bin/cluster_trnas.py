#!/usr/bin/env python3
"""Replace clustering.pl: cluster identical mature tRNA sequences.

Usage: cluster_trnas.py <input.fa> <cluster.fa> <clusterInfo.fa>

Logic (matching legacy clustering.pl):
  - Read single-line FASTA (one header + one sequence per entry)
  - Group entries with identical sequences
  - Output cluster.fa: one entry per unique sequence, named >clusterN
  - Output clusterInfo.fa: all original entries annotated with their cluster ID
    as >clusterN:original_header_without_>

Note: the legacy Perl script iterates hash keys in non-deterministic order,
so cluster numbering varies between runs. This Python version uses insertion
order (first-seen sequence gets cluster1), making output deterministic.
"""

import sys
from collections import OrderedDict


def cluster_trnas(in_path, cluster_path, info_path):
    # Read single-line FASTA
    entries = []
    with open(in_path) as fh:
        lines = [l.rstrip("\n") for l in fh]

    for i in range(0, len(lines), 2):
        header = lines[i]
        seq = lines[i + 1]
        entries.append((header, seq))

    # Group by sequence (preserving first-seen order)
    seq_to_headers = OrderedDict()
    for header, seq in entries:
        if seq not in seq_to_headers:
            seq_to_headers[seq] = []
        seq_to_headers[seq].append(header)

    # Write outputs
    with open(cluster_path, "w") as cluster_fh, open(info_path, "w") as info_fh:
        for cluster_no, (seq, headers) in enumerate(seq_to_headers.items(), start=1):
            cluster_fh.write(f">cluster{cluster_no}\n{seq}\n")
            for header in headers:
                # Strip leading '>' from original header
                bare_header = header.lstrip(">")
                info_fh.write(f">cluster{cluster_no}:{bare_header}\n{seq}\n")


if __name__ == "__main__":
    cluster_trnas(sys.argv[1], sys.argv[2], sys.argv[3])
