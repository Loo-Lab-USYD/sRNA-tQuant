#!/usr/bin/env python3
"""Convert Salmon quant.sf to anticodon-level RPM with optional span-weighted blending.

Without --bam: behaves identically to salmon_to_rpm.py (backward compatible).

With --bam: blends Salmon EM estimates with direct counts from reads spanning
the anticodon position (nucleotides 34-36), improving isoacceptor resolution
for amino acid families with near-identical tRNA sequences.

Usage:
    # Standard mode (same as salmon_to_rpm.py):
    span_weighted_rpm.py quant.sf clusterInfo.fa output.csv

    # Span-weighted mode:
    span_weighted_rpm.py quant.sf clusterInfo.fa output.csv \\
        --bam topscore.bam \\
        --annotations sample.annotations.csv

    # With all options:
    span_weighted_rpm.py quant.sf clusterInfo.fa output.csv \\
        --bam topscore.bam \\
        --cluster-rpm sample.cluster_rpm.csv \\
        --annotations sample.annotations.csv \\
        --min-spanning 100 \\
        --max-blend-weight 0.8
"""

import argparse
import re
import sys
from collections import defaultdict

import pandas as pd


# ---------------------------------------------------------------------------
# Reference parsing
# ---------------------------------------------------------------------------

def parse_cluster_info(info_path):
    """Parse clusterInfo.fa -> {cluster: anticodon} mapping + sequences."""
    mapping = {}
    sequences = {}
    ac_pattern = re.compile(r"(?<=-)i?[A-Z]{1}[a-zC]+[ACTGN]{3}")
    current_id = None

    with open(info_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cid = line[1:].split(":")[0]
                m = ac_pattern.search(line)
                if m:
                    mapping[cid] = m.group()
                current_id = cid
                sequences[cid] = ""
            elif current_id:
                sequences[current_id] += line

    return mapping, sequences


def get_amino(ac):
    """Extract amino acid from anticodon name: ThrAGT -> Thr, iMetCAT -> iMet."""
    m = re.match(r"(i?[A-Z][a-z]+)", ac)
    return m.group(1) if m else ac


def compute_sibling_similarity(mapping, sequences):
    """For each anticodon, compute max sequence identity to sibling isoacceptors.

    Siblings = anticodons of the same amino acid. Returns {anticodon: float}.
    """
    # Group clusters by amino acid -> anticodon
    amino_groups = defaultdict(lambda: defaultdict(list))
    for cid, ac in mapping.items():
        amino_groups[get_amino(ac)][ac].append(cid)

    sibling_sim = {}
    for amino, ac_clusters in amino_groups.items():
        anticodons = list(ac_clusters.keys())
        for ac1 in anticodons:
            best = 0.0
            for ac2 in anticodons:
                if ac1 == ac2:
                    continue
                for c1 in ac_clusters[ac1]:
                    for c2 in ac_clusters[ac2]:
                        s1 = sequences.get(c1, "")
                        s2 = sequences.get(c2, "")
                        if s1 and s2:
                            maxlen = max(len(s1), len(s2))
                            matches = sum(a == b for a, b in zip(s1, s2))
                            best = max(best, matches / maxlen if maxlen > 0 else 0)
            sibling_sim[ac1] = best

    return sibling_sim


# ---------------------------------------------------------------------------
# BAM spanning-read analysis
# ---------------------------------------------------------------------------

AC_START = 34  # 1-indexed, first anticodon nucleotide
AC_END = 36    # 1-indexed, third anticodon nucleotide


def _parse_cigar_ref_end(cigar, pos):
    """Compute reference end position from CIGAR string and 1-based start."""
    ref_pos = pos
    for length_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in ("M", "D", "N", "=", "X"):
            ref_pos += int(length_str)
    return ref_pos - 1  # inclusive end


def count_spanning_reads(bam_path, mapping):
    """Count reads spanning the anticodon (pos 34-36) per anticodon.

    Uses pysam if available, falls back to samtools subprocess.
    Returns (ac_spanning, ac_total) dicts: anticodon -> count.
    """
    try:
        import pysam
        return _count_spanning_pysam(bam_path, mapping)
    except ImportError:
        return _count_spanning_samtools(bam_path, mapping)


def _count_spanning_pysam(bam_path, mapping):
    import pysam

    ac_spanning = defaultdict(int)
    ac_total = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            ref_name = read.reference_name
            if ref_name not in mapping:
                continue

            ac = mapping[ref_name]
            ac_total[ac] += 1

            ref_start = read.reference_start + 1  # 0-based -> 1-based
            ref_end = read.reference_end  # pysam exclusive end = 1-based inclusive

            if ref_start <= AC_START and ref_end >= AC_END:
                ac_spanning[ac] += 1

    return dict(ac_spanning), dict(ac_total)


def _count_spanning_samtools(bam_path, mapping):
    import subprocess

    ac_spanning = defaultdict(int)
    ac_total = defaultdict(int)

    proc = subprocess.Popen(
        ["samtools", "view", bam_path],
        stdout=subprocess.PIPE, text=True,
    )
    for line in proc.stdout:
        fields = line.split("\t", 6)
        ref_name = fields[2]
        if ref_name not in mapping:
            continue

        ac = mapping[ref_name]
        ac_total[ac] += 1

        pos = int(fields[3])  # 1-based
        cigar = fields[5]
        ref_end = _parse_cigar_ref_end(cigar, pos)

        if pos <= AC_START and ref_end >= AC_END:
            ac_spanning[ac] += 1

    proc.wait()
    return dict(ac_spanning), dict(ac_total)


# ---------------------------------------------------------------------------
# Blending logic
# ---------------------------------------------------------------------------

def compute_blend_weights(anticodons, sibling_sim, ac_spanning, ac_total,
                          min_spanning=100, max_weight=0.8):
    """Compute per-anticodon blend weight (0 = all Salmon, 1 = all spanning).

    Rules:
        - Single-isoacceptor amino acids: w = 0 (Salmon is unambiguous)
        - Sibling similarity < 0.85: w = 0 (isoacceptors distinguishable)
        - Any isoacceptor in the family has < min_spanning reads: w = 0
        - Otherwise: w = sibling_sim × (1 - resolvability), capped at max_weight
    """
    # Group anticodons by amino acid
    amino_groups = defaultdict(list)
    for ac in anticodons:
        amino_groups[get_amino(ac)].append(ac)

    weights = {}
    for ac in anticodons:
        aa = get_amino(ac)
        siblings = amino_groups[aa]

        # Single-isoacceptor: no blending needed
        if len(siblings) <= 1:
            weights[ac] = 0.0
            continue

        sim = sibling_sim.get(ac, 0.0)

        # Low similarity: Salmon can distinguish
        if sim < 0.85:
            weights[ac] = 0.0
            continue

        # Check spanning read depth for all isoacceptors in this family
        family_has_enough = all(
            ac_spanning.get(sib, 0) >= min_spanning for sib in siblings
        )
        if not family_has_enough:
            weights[ac] = 0.0
            continue

        # Compute resolvability for this anticodon
        total = ac_total.get(ac, 0)
        spanning = ac_spanning.get(ac, 0)
        resolvability = spanning / total if total > 0 else 0

        w = sim * (1 - resolvability)
        weights[ac] = min(w, max_weight)

    return weights


# ---------------------------------------------------------------------------
# Main quantification
# ---------------------------------------------------------------------------

def salmon_to_rpm(quant_path, info_path, out_path, cluster_rpm_path=None,
                  bam_path=None, annotations_path=None,
                  min_spanning=100, max_blend_weight=0.8, clamp_threshold=0):
    """Compute anticodon-level RPMs with optional span-weighted blending."""

    # Parse reference
    mapping, sequences = parse_cluster_info(info_path)
    anticodons = sorted(set(mapping.values()))

    # Load Salmon quantification
    quant = pd.read_table(quant_path)
    tot_reads = quant["NumReads"].sum()

    if tot_reads == 0:
        rpm_df = pd.DataFrame(0.0, index=anticodons, columns=["rpm", "rpm_nomod"])
        reads_row = pd.DataFrame({"rpm": [0], "rpm_nomod": [0]}, index=["READS"])
        rpm_df = pd.concat([rpm_df, reads_row])
        rpm_df.to_csv(out_path)
        if cluster_rpm_path:
            pd.DataFrame(columns=["cluster", "rpm"]).to_csv(cluster_rpm_path, index=False)
        return

    # Salmon per-anticodon RPM
    salmon_rpm = pd.Series(0.0, index=anticodons)
    cluster_rpms = {}
    for _, row in quant.iterrows():
        cluster_name = row["Name"]
        if cluster_name not in mapping:
            continue
        ac = mapping[cluster_name]
        crpm = row["NumReads"] * 1_000_000 / tot_reads
        salmon_rpm[ac] += crpm
        cluster_rpms[cluster_name] = crpm

    # Without BAM: standard mode (identical to salmon_to_rpm.py)
    if not bam_path:
        rpm_df = pd.DataFrame({"rpm": salmon_rpm, "rpm_nomod": salmon_rpm})
        # Clamp low-abundance anticodons to zero (removes tRNA fragment noise)
        if clamp_threshold > 0:
            for col in ["rpm", "rpm_nomod"]:
                rpm_df.loc[rpm_df[col] < clamp_threshold, col] = 0.0
        reads_row = pd.DataFrame(
            {"rpm": [tot_reads], "rpm_nomod": [tot_reads]}, index=["READS"]
        )
        rpm_df = pd.concat([rpm_df, reads_row])
        rpm_df.to_csv(out_path)

        if cluster_rpm_path:
            cluster_df = pd.DataFrame(
                [{"cluster": k, "rpm": v} for k, v in cluster_rpms.items()]
            )
            cluster_df.to_csv(cluster_rpm_path, index=False)
        return

    # --- Span-weighted mode ---

    # Count spanning reads
    print(f"  Counting anticodon-spanning reads from {bam_path}...", file=sys.stderr)
    ac_spanning, ac_total = count_spanning_reads(bam_path, mapping)

    total_spanning = sum(ac_spanning.values())
    print(f"  Spanning reads: {total_spanning:,} / {sum(ac_total.values()):,} "
          f"({100 * total_spanning / max(sum(ac_total.values()), 1):.1f}%)",
          file=sys.stderr)

    # Spanning-only RPM
    span_rpm = pd.Series(0.0, index=anticodons)
    if total_spanning > 0:
        for ac in anticodons:
            span_rpm[ac] = ac_spanning.get(ac, 0) * 1_000_000 / total_spanning

    # Compute sibling similarity (reference-level, deterministic)
    sibling_sim = compute_sibling_similarity(mapping, sequences)

    # Compute blend weights
    weights = compute_blend_weights(
        anticodons, sibling_sim, ac_spanning, ac_total,
        min_spanning=min_spanning, max_weight=max_blend_weight,
    )

    # Blend
    blended_rpm = pd.Series(0.0, index=anticodons)
    for ac in anticodons:
        w = weights[ac]
        blended_rpm[ac] = w * span_rpm[ac] + (1 - w) * salmon_rpm[ac]

    # Output RPM (blended values in both columns for compatibility)
    rpm_df = pd.DataFrame({"rpm": blended_rpm, "rpm_nomod": blended_rpm})
    # Clamp low-abundance anticodons to zero (removes tRNA fragment noise)
    if clamp_threshold > 0:
        for col in ["rpm", "rpm_nomod"]:
            rpm_df.loc[rpm_df[col] < clamp_threshold, col] = 0.0
    reads_row = pd.DataFrame(
        {"rpm": [tot_reads], "rpm_nomod": [tot_reads]}, index=["READS"]
    )
    rpm_df = pd.concat([rpm_df, reads_row])
    rpm_df.to_csv(out_path)

    # Cluster RPM (unchanged — always from Salmon)
    if cluster_rpm_path:
        cluster_df = pd.DataFrame(
            [{"cluster": k, "rpm": v} for k, v in cluster_rpms.items()]
        )
        cluster_df.to_csv(cluster_rpm_path, index=False)

    # Annotations
    if annotations_path:
        ann_rows = []
        for ac in anticodons:
            total = ac_total.get(ac, 0)
            spanning = ac_spanning.get(ac, 0)
            resolvability = spanning / total if total > 0 else 0
            sim = sibling_sim.get(ac, 0.0)
            w = weights[ac]
            n_siblings = sum(1 for a in anticodons if get_amino(a) == get_amino(ac)) - 1

            if n_siblings == 0:
                tier = "high"
            elif sim < 0.85:
                tier = "high"
            elif w > 0:
                tier = "moderate"
            else:
                tier = "low"

            ann_rows.append({
                "anticodon": ac,
                "salmon_rpm": salmon_rpm[ac],
                "spanning_rpm": span_rpm[ac],
                "blended_rpm": blended_rpm[ac],
                "blend_weight": w,
                "n_spanning_reads": spanning,
                "n_total_reads": total,
                "resolvability": resolvability,
                "sibling_similarity": sim,
                "n_sibling_anticodons": n_siblings,
                "confidence_tier": tier,
            })
        ann_df = pd.DataFrame(ann_rows).set_index("anticodon")
        ann_df.to_csv(annotations_path)
        n_blended = sum(1 for w in weights.values() if w > 0)
        print(f"  Blended {n_blended}/{len(anticodons)} anticodons "
              f"(annotations: {annotations_path})", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Salmon quant.sf to anticodon RPM with optional span-weighted blending"
    )
    parser.add_argument("quant", help="Salmon quant.sf file")
    parser.add_argument("info", help="clusterInfo.fa file")
    parser.add_argument("output", help="Output RPM CSV")
    parser.add_argument("--bam", help="Topscore BAM for span-weighted mode")
    parser.add_argument("--cluster-rpm", help="Output per-cluster RPM CSV")
    parser.add_argument("--annotations", help="Output per-anticodon annotation CSV")
    parser.add_argument("--min-spanning", type=int, default=5,
                        help="Minimum spanning reads per isoacceptor to enable blending (default: 5)")
    parser.add_argument("--max-blend-weight", type=float, default=0.55,
                        help="Maximum blend weight for spanning reads (default: 0.55)")
    parser.add_argument("--clamp-threshold", type=float, default=0,
                        help="Zero out anticodons below this RPM (default: 0 = disabled)")
    args = parser.parse_args()

    salmon_to_rpm(
        args.quant, args.info, args.output,
        cluster_rpm_path=getattr(args, "cluster_rpm"),
        bam_path=args.bam,
        annotations_path=args.annotations,
        min_spanning=args.min_spanning,
        max_blend_weight=args.max_blend_weight,
        clamp_threshold=args.clamp_threshold,
    )
