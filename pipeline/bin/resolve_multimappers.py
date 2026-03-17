#!/usr/bin/env python3
"""Replace distribute_reads.py: anticodon-aware multimapper resolution.

Usage: resolve_multimappers.py <input.sam> <clusterInfo.fa> <output.sam> <multimappers.sam>

Logic (matching legacy distribute_reads.py):
  - Parse clusterInfo.fa to build cluster_id → anticodon mapping
  - Stream name-sorted SAM:
    - SAM headers (@) → write to output
    - Unique reads (NH:i:1) → write to output
    - Multimappers: group consecutive lines by read name
      - If all hits share the same anticodon → write one alignment to output
      - If hits span different anticodons → write all to multimappers file

Note: the legacy code has a bug where the final group of multimappers is never
flushed (loop ends before the "new read name" branch processes the last group).
This version fixes that by flushing after the loop.
"""

import re
import sys


def build_cluster_mapping(info_path):
    """Parse clusterInfo.fa headers to build cluster → anticodon mapping."""
    mapping = {}
    with open(info_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cluster_id = re.findall(r"cluster[0-9]+", line)[0]
                # Extract anticodon identifier: e.g. AlaAGC, iMetCAT, SerCGA
                anticodon = re.findall(r"(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}", line)[0]
                mapping[cluster_id] = anticodon
    return mapping


def resolve_multimappers(sam_in, info_path, sam_out, multi_out):
    mapping = build_cluster_mapping(info_path)

    with open(sam_in) as in_fh, \
         open(sam_out, "w") as out_fh, \
         open(multi_out, "w") as multi_fh:

        # State for grouping multimappers
        current_read = None
        anticodons = []
        lines = []

        def flush_group():
            """Process accumulated multimapper group."""
            if not lines:
                return
            if len(set(anticodons)) == 1:
                # All hits share same anticodon → write last alignment
                out_fh.write(lines[-1])
            elif len(set(anticodons)) > 1:
                # Different anticodons → write all to multimappers
                for item in lines:
                    multi_fh.write(item)

        for line in in_fh:
            if line[0] == "@":
                out_fh.write(line)
            elif "NH:i:1\t" in line:
                out_fh.write(line)
            else:
                fields = line.split("\t")
                read_name = fields[0]
                cluster_id = fields[2]

                if read_name == current_read:
                    anticodons.append(mapping[cluster_id])
                    lines.append(line)
                else:
                    flush_group()
                    current_read = read_name
                    anticodons = [mapping[cluster_id]]
                    lines = [line]

        # Flush final group (fixes legacy bug)
        flush_group()


if __name__ == "__main__":
    resolve_multimappers(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
