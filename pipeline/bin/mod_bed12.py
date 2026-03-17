#!/usr/bin/env python3
"""Add flanking regions to tRNA BED12 coordinates for pre-tRNA extraction.

Usage: mod_bed12.py <input.bed12> <output.bed12> [--flank N]

Replaces modBed12.pl. Extends BED12 start/end by the flanking amount
(default 50nt) to capture pre-tRNA sequences including leader/trailer.
Block sizes and starts are updated so that bedtools getfasta -split
includes the flanking regions in the extracted sequences.

Pseudogenes (name containing 'pseudo') are excluded (matching legacy).
"""

import argparse
import sys


def mod_bed12(in_path, out_path, flank=50):
    with open(in_path) as fh_in, open(out_path, "w") as fh_out:
        for line in fh_in:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                continue

            # Skip pseudogenes (matching legacy modBed12.pl)
            if "pseudo" in fields[3].lower():
                continue

            chrom = fields[0]
            orig_start = int(fields[1])
            orig_end = int(fields[2])
            name = fields[3]
            score = fields[4]
            strand = fields[5]
            rgb = fields[8]
            block_count = int(fields[9])
            block_sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
            block_starts = [int(x) for x in fields[11].rstrip(",").split(",")]

            new_start = max(0, orig_start - flank)
            actual_upstream = orig_start - new_start  # may be < flank if near chrom start
            new_end = orig_end + flank
            orig_span = orig_end - orig_start

            if block_count == 2:
                # Spliced tRNA (has intron): extend first block upstream,
                # last block downstream, recalculate second block start
                # Legacy: blockSizes = [s1+50, s2+50], blockStarts = [0, span-s2+50]
                new_block_sizes = [
                    block_sizes[0] + actual_upstream,
                    block_sizes[1] + flank,
                ]
                new_block_starts = [
                    0,
                    (orig_span - block_sizes[1]) + actual_upstream,
                ]
            else:
                # Single-block tRNA: extend block size by both flanks
                # Legacy: blockSize = original + 100 (50 each side)
                new_block_sizes = [block_sizes[0] + actual_upstream + flank]
                new_block_starts = [0]

            sizes_str = ",".join(str(s) for s in new_block_sizes) + ","
            starts_str = ",".join(str(s) for s in new_block_starts) + ","

            fh_out.write(
                f"{chrom}\t{new_start}\t{new_end}\t{name}\t{score}\t{strand}\t"
                f"{new_start}\t{new_end}\t{rgb}\t{block_count}\t"
                f"{sizes_str}\t{starts_str}\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add flanking regions to BED12")
    parser.add_argument("input", help="Input BED12 file")
    parser.add_argument("output", help="Output BED12 file")
    parser.add_argument("--flank", type=int, default=50, help="Flanking size in nt (default: 50)")
    args = parser.parse_args()
    mod_bed12(args.input, args.output, args.flank)
