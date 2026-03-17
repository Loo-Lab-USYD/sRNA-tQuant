#!/usr/bin/env python3
"""Replace extract_data.py + extract_modifications.py: merge per-sample results.

Usage: aggregate_results.py <output_dir> <rpm_csv1> [rpm_csv2 ...] [--ase ase_csv1 ase_csv2 ...]

Simplified replacement for the legacy scripts which had complex TCGA-specific
sample naming logic. In the Nextflow pipeline, sample IDs come from the
samplesheet and RPM files are named consistently.

Outputs:
  <output_dir>/results_rpm.csv        - wide matrix (anticodons × samples) with rpm values
  <output_dir>/results_rpm_nomod.csv  - wide matrix with rpm_nomod values
  <output_dir>/modifications.csv      - merged modification calls (if ASE files provided)
"""

import os
import sys

import pandas as pd


def aggregate_results(output_dir, rpm_files, ase_files=None):
    os.makedirs(output_dir, exist_ok=True)

    # Merge RPM files
    mod_frames = {}
    nomod_frames = {}

    for rpm_file in rpm_files:
        # Derive sample ID from filename: /path/to/sampleA.RPM.csv → sampleA
        sample_id = os.path.basename(rpm_file).replace(".RPM.csv", "")
        df = pd.read_csv(rpm_file, index_col=0)
        mod_frames[sample_id] = df["rpm"]
        nomod_frames[sample_id] = df["rpm_nomod"]

    results_mod = pd.DataFrame(mod_frames)
    results_nomod = pd.DataFrame(nomod_frames)

    results_mod.to_csv(os.path.join(output_dir, "results_rpm.csv"))
    results_nomod.to_csv(os.path.join(output_dir, "results_rpm_nomod.csv"))

    # Merge ASE/modification files if provided
    if ase_files:
        merged = []
        for ase_file in ase_files:
            sample_id = os.path.basename(ase_file).replace(".ASE.csv", "")
            df = pd.read_csv(ase_file, sep="\t")
            df["sample"] = sample_id
            df["modif_id"] = df.apply(
                lambda x: f"{x['contig']}-{x['position']}-{x['refAllele']}to{x['altAllele']}",
                axis=1,
            )
            merged.append(df)

        if merged:
            modifications = pd.concat(merged, ignore_index=True)
            modifications.to_csv(
                os.path.join(output_dir, "modifications.csv"), index=False
            )


if __name__ == "__main__":
    args = sys.argv[1:]
    output_dir = args[0]

    # Split args on --ase flag
    if "--ase" in args:
        ase_idx = args.index("--ase")
        rpm_files = args[1:ase_idx]
        ase_files = args[ase_idx + 1 :]
    else:
        rpm_files = args[1:]
        ase_files = None

    aggregate_results(output_dir, rpm_files, ase_files)
