#!/usr/bin/env python3
"""Adjust RPM values based on post-transcriptional modifications in the anticodon.

Usage: adjust_modifications.py <rpm.csv> <ase.csv> <clusterInfo.fa> <tRNAs.txt> <output.csv> [--cluster-rpm <cluster_rpm.csv>]

Logic:
  - Takes the base RPM CSV (from Salmon or Decision quantification)
  - Loads ASE (allele-specific expression) data from variant calling
  - For each cluster with variants in the anticodon triplet:
    redistribute reads proportionally based on alt/ref allele ratios
  - When --cluster-rpm is provided, uses per-cluster RPMs for accurate
    redistribution (fixes bug where total anticodon RPM was used as base)
  - Outputs adjusted RPM CSV
"""

import argparse
import re

import pandas as pd


def adjust_modifications(rpm_path, ase_path, info_path, trnas_path, out_path, cluster_rpm_path=None):
    # Load base RPM
    rpm = pd.read_csv(rpm_path, index_col=0)

    # Load ASE data
    ase = pd.read_table(ase_path, index_col="contig")

    if ase.empty:
        # No variants found — output unchanged
        rpm.to_csv(out_path)
        return

    # Load cluster info for sequence and anticodon mapping
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
        trnas = [t for t in fh.read().split("\n") if t]

    # Load per-cluster RPMs if available
    cluster_rpms = {}
    if cluster_rpm_path:
        cdf = pd.read_csv(cluster_rpm_path)
        cluster_rpms = dict(zip(cdf["cluster"], cdf["rpm"]))

    # Reset rpm column (keep rpm_nomod as-is)
    for trna in trnas:
        if trna in rpm.index:
            rpm.loc[trna, "rpm"] = rpm.loc[trna, "rpm_nomod"]

    # Process each cluster that has ASE variants
    for cluster in set(ase.index):
        if cluster not in mapping:
            continue

        anticodon_label = mapping[cluster]
        if anticodon_label == "UndetNNN":
            continue

        # Localize anticodon triplet in the cluster sequence
        anticodon_3 = anticodon_label[-3:]
        seq = cluster_seq.get(cluster, "")
        if not seq:
            continue

        acod = re.findall(anticodon_3, seq, re.I)
        match = re.search(anticodon_3, seq, re.I)
        if not match:
            continue
        pos = match.start() + 1

        # Disambiguate if anticodon triplet appears multiple times
        trim = 0
        end5 = 0
        end3 = 1
        while len(acod) != 1:
            substr = seq[(trim + end5): -(trim + end3)]
            acod = re.findall(anticodon_3, substr, re.I)
            m = re.search(anticodon_3, substr, re.I)
            if not m:
                break
            pos = m.start() + trim + 1
            if end3 == 0:
                end3 += 1
            elif end5 == 0:
                end5 += 1
            else:
                trim += 1
                end5 = 0
                end3 = 0

        # Check if ASE variants fall within the anticodon
        ase_c = ase.loc[[cluster], :]
        ase_acod = ase_c.loc[
            [p in range(pos, pos + 3) for p in ase_c.position],
        ]

        if ase_acod.empty:
            continue

        # Calculate the fraction of reads with the alternative allele
        total_depth = sum(ase_acod.loc[:, "rawDepth"])
        if total_depth == 0:
            continue
        alt_fraction = sum(ase_acod.loc[:, "altCount"]) / total_depth

        # Get the RPM for this cluster (per-cluster if available, else total anticodon)
        if cluster in cluster_rpms:
            base_rpm = cluster_rpms[cluster]
        else:
            base_rpm = rpm.loc[anticodon_label, "rpm_nomod"]
        alt_rpms = base_rpm * alt_fraction

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

        # Redistribute
        if new_trna in rpm.index:
            rpm.loc[new_trna, "rpm"] += alt_rpms
        rpm.loc[anticodon_label, "rpm"] -= alt_rpms

    rpm.to_csv(out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adjust RPM for anticodon modifications")
    parser.add_argument("rpm", help="Input RPM CSV")
    parser.add_argument("ase", help="ASE CSV from count_alleles")
    parser.add_argument("info", help="clusterInfo.fa")
    parser.add_argument("trnas", help="tRNAs.txt")
    parser.add_argument("output", help="Output adjusted RPM CSV")
    parser.add_argument("--cluster-rpm", help="Per-cluster RPM CSV from salmon_to_rpm")
    args = parser.parse_args()

    adjust_modifications(args.rpm, args.ase, args.info, args.trnas, args.output,
                         getattr(args, "cluster_rpm"))
