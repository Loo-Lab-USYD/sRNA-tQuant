#!/usr/bin/env python3
"""Convert Salmon quant.sf to anticodon-level RPM CSV.

Usage: salmon_to_rpm.py <quant.sf> <clusterInfo.fa> <output.csv> [--cluster-rpm <cluster_rpm.csv>]

Logic:
  - Parse Salmon's quant.sf (columns: Name, Length, EffectiveLength, TPM, NumReads)
  - Map each cluster to its anticodon via clusterInfo.fa headers
  - Sum NumReads per anticodon
  - Convert to RPM (reads per million mapped reads)
  - Output CSV with columns: [anticodon, rpm, rpm_nomod]
    (rpm_nomod == rpm here; modification adjustment is a separate step)
  - Optionally output per-cluster RPMs for downstream modification adjustment
"""

import argparse
import re

import pandas as pd


def salmon_to_rpm(quant_path, info_path, out_path, cluster_rpm_path=None):
    # Load Salmon quantification
    quant = pd.read_table(quant_path)

    # Load cluster→anticodon mapping from clusterInfo headers
    with open(info_path) as fh:
        raw = fh.read().split("\n")
    headers = raw[0:-1:2]

    mapping = {}
    for s in headers:
        cluster = re.findall(r"cluster[0-9]+", s)[0]
        anticodon = re.findall(r"(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}", s)[0]
        mapping[cluster] = anticodon

    # Get unique anticodons for the index
    anticodons = sorted(set(mapping.values()))

    # Sum NumReads per anticodon
    rpm_df = pd.DataFrame(0.0, index=anticodons, columns=["rpm", "rpm_nomod"])
    tot_reads = quant["NumReads"].sum()

    if tot_reads == 0:
        # No reads mapped — write empty RPM
        reads_row = pd.DataFrame(
            {"rpm": [0], "rpm_nomod": [0]}, index=["READS"]
        )
        rpm_df = pd.concat([rpm_df, reads_row])
        rpm_df.to_csv(out_path)
        if cluster_rpm_path:
            pd.DataFrame(columns=["cluster", "rpm"]).to_csv(cluster_rpm_path, index=False)
        return

    cluster_rpms = {}
    for _, row in quant.iterrows():
        cluster_name = row["Name"]
        if cluster_name not in mapping:
            continue
        anticodon = mapping[cluster_name]
        cluster_rpm = row["NumReads"] * 1_000_000 / tot_reads
        rpm_df.loc[anticodon, "rpm"] += cluster_rpm
        rpm_df.loc[anticodon, "rpm_nomod"] += cluster_rpm
        cluster_rpms[cluster_name] = cluster_rpm

    # Append total reads row
    reads_row = pd.DataFrame(
        {"rpm": [tot_reads], "rpm_nomod": [tot_reads]}, index=["READS"]
    )
    rpm_df = pd.concat([rpm_df, reads_row])

    rpm_df.to_csv(out_path)

    # Output per-cluster RPMs if requested
    if cluster_rpm_path:
        cluster_df = pd.DataFrame(
            [{"cluster": k, "rpm": v} for k, v in cluster_rpms.items()]
        )
        cluster_df.to_csv(cluster_rpm_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Salmon quant.sf to anticodon-level RPM")
    parser.add_argument("quant", help="Salmon quant.sf file")
    parser.add_argument("info", help="clusterInfo.fa file")
    parser.add_argument("output", help="Output RPM CSV")
    parser.add_argument("--cluster-rpm", help="Output per-cluster RPM CSV")
    args = parser.parse_args()

    salmon_to_rpm(args.quant, args.info, args.output, getattr(args, "cluster_rpm"))
