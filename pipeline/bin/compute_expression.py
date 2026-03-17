#!/usr/bin/env python3
"""Replace compute_expression.py: calculate RPM per anticodon with optional modification adjustment.

Usage: compute_expression.py <expression.txt> <ASE.csv> <clusterInfo.fa> <tRNAs.txt> <output.csv>

Logic (matching legacy compute_expression.py):
  - Load idxstats expression counts (cluster → read count)
  - Load ASE (allele-specific expression) variant calls
  - Map clusters to anticodons via clusterInfo headers
  - For each cluster:
    - Always add RPM to rpm_nomod column
    - If no ASE variants, or variants are outside the anticodon triplet:
      add RPM to rpm column unchanged
    - If ASE variants fall within the anticodon triplet:
      redistribute a fraction of reads to the alternative anticodon
      based on altCount/rawDepth ratio
  - Append a READS row with total mapped read count
  - Output CSV with columns: [anticodon, rpm, rpm_nomod]

Fix vs legacy: the READS row index is dynamically determined instead of
hardcoded to position 66.
"""

import re
import sys

import pandas as pd


def compute_expression(expr_path, ase_path, info_path, trnas_path, out_path):
    # Load expression data (idxstats format: contig, length, mapped, unmapped)
    expression = pd.read_table(expr_path, header=None, index_col=0)
    expression.drop("*", inplace=True)
    expression = expression[expression.loc[:, 2] != 0]

    # Load ASE data
    ase = pd.read_table(ase_path, index_col="contig")

    # Load cluster info
    with open(info_path) as fh:
        raw = fh.read().split("\n")
    seqs = raw[1::2]
    headers = raw[0:-1:2]

    cluster_seq = dict(
        zip([re.findall(r"cluster[0-9]+", s)[0] for s in headers], seqs)
    )
    mapping = {
        re.findall(r"cluster[0-9]+", s)[0]: re.findall(
            r"(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}", s
        )[0]
        for s in headers
    }

    # Load tRNA list
    with open(trnas_path) as fh:
        trnas = fh.read().split("\n")
    # Remove empty trailing entry if present
    trnas = [t for t in trnas if t]

    # Initialize RPM dataframe
    rpm = pd.DataFrame(0.0, index=trnas, columns=["rpm", "rpm_nomod"])
    tot_reads = sum(expression.loc[:, 2])

    # Compute expression per cluster
    for c in expression.index:
        cluster_rpm = expression.loc[c, 2] * 1_000_000 / tot_reads
        rpm.loc[mapping[c], "rpm_nomod"] += cluster_rpm

        if c not in ase.index:
            rpm.loc[mapping[c], "rpm"] += cluster_rpm
        elif mapping[c] == "UndetNNN":
            rpm.loc[mapping[c], "rpm"] += cluster_rpm
        else:
            # Localize anticodon triplet in the cluster sequence
            anticodon_3 = mapping[c][-3:]
            seq = cluster_seq[c]
            acod = re.findall(anticodon_3, seq, re.I)
            pos = re.search(anticodon_3, seq, re.I).start() + 1

            # If anticodon triplet appears multiple times, trim from ends
            trim = 0
            end5 = 0
            end3 = 1
            while len(acod) != 1:
                substr = seq[(trim + end5) : -(trim + end3)]
                acod = re.findall(anticodon_3, substr, re.I)
                pos = re.search(anticodon_3, substr, re.I).start() + trim + 1
                if end3 == 0:
                    end3 += 1
                elif end5 == 0:
                    end5 += 1
                else:
                    trim += 1
                    end5 = 0
                    end3 = 0

            # Check if ASE variants fall within the anticodon
            ase_c = ase.loc[[c], :]
            ase_acod = ase_c.loc[
                [p in range(pos, pos + 3) for p in ase_c.position],
            ]

            if ase_acod.empty:
                # Modification does not affect anticodon
                rpm.loc[mapping[c], "rpm"] += cluster_rpm
            else:
                # Redistribute reads based on modification frequency
                alt_rpms = cluster_rpm * (
                    sum(ase_acod.loc[:, "altCount"])
                    / sum(ase_acod.loc[:, "rawDepth"])
                )

                # Determine the alternative anticodon
                new_acod = list(acod[0])
                for p in ase_acod.position:
                    acod_pos = p - pos
                    new_acod[acod_pos] = ase_acod.loc[
                        ase_acod.position == p, "altAllele"
                    ].iloc[0]
                new_acod = "".join(new_acod).upper()

                new_trna = [s for s in trnas if new_acod in s]
                new_trna = new_trna[0] if new_trna else "UndetNNN"

                rpm.loc[new_trna, "rpm"] += alt_rpms
                rpm.loc[mapping[c], "rpm"] += cluster_rpm - alt_rpms

    # Append total reads row (dynamic, not hardcoded to index 66)
    reads_row = pd.DataFrame(
        {"rpm": [tot_reads], "rpm_nomod": [tot_reads]}, index=["READS"]
    )
    rpm = pd.concat([rpm, reads_row])

    rpm.to_csv(out_path)


if __name__ == "__main__":
    compute_expression(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
