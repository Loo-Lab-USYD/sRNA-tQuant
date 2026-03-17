#!/usr/bin/env python3
"""
Validate K562 ENCODE integration test outputs from the tRNA mapping pipeline.

Extends the base validator with K562-specific checks:
  - Correct sample IDs (K562_rep1_test, K562_rep2_test)
  - Replicate correlation (both samples should be very similar)
  - fastp handled 101bp → small RNA trimming correctly

Usage:
    python3 validate_k562_outputs.py <results_dir>
"""

import csv
import json
import re
import sys
from pathlib import Path


ANTICODON_PATTERN = re.compile(r'^i?[A-Z][a-zA-Z()]+[ACTGN]{3}$')

K562_SAMPLES = ['K562_rep1_test', 'K562_rep2_test']


def validate(results_dir: Path, work_dir: Path | None = None) -> list[str]:
    """Run all validations, return list of failure messages."""
    failures = []

    # --- results_rpm.csv ---
    rpm_path = results_dir / 'results_rpm.csv'
    if not rpm_path.exists():
        failures.append('results_rpm.csv does not exist')
        return failures

    with open(rpm_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    if len(header) < 3:
        failures.append(f'RPM matrix has {len(header)} columns, expected >= 3')

    if len(rows) < 3:
        failures.append(f'RPM matrix has {len(rows)} data rows, expected >= 3')

    # Check K562 sample IDs in header
    sample_ids = header[1:]
    for expected in K562_SAMPLES:
        if not any(expected in s for s in sample_ids):
            failures.append(f'Sample {expected} not in RPM header: {sample_ids}')

    # Check anticodon names
    index_values = [row[0] for row in rows]
    anticodon_rows = [v for v in index_values if v != 'READS']

    if 'READS' not in index_values:
        failures.append('READS row not found in RPM matrix')

    bad_names = []
    for name in anticodon_rows:
        if not ANTICODON_PATTERN.match(name):
            bad_names.append(name)
    if bad_names:
        failures.append(f'Invalid anticodon names: {bad_names[:5]}')

    # Check no raw FASTA headers
    for name in index_values:
        if '::' in name or '>' in name:
            failures.append(f'Raw FASTA header in index: {name}')
            break

    # Check RPM values are non-negative
    for row in rows:
        for val in row[1:]:
            try:
                v = float(val)
                if v < 0:
                    failures.append(f'Negative RPM value: {val} in row {row[0]}')
                    break
            except ValueError:
                failures.append(f'Non-numeric RPM value: {val!r} in row {row[0]}')
                break

    # Check READS row has positive values
    reads_rows = [row for row in rows if row[0] == 'READS']
    if reads_rows:
        reads_vals = [float(v) for v in reads_rows[0][1:]]
        if not all(v > 0 for v in reads_vals):
            failures.append(f'READS row has zero/negative values: {reads_vals}')
        min_reads = min(reads_vals)
        if min_reads < 10:
            failures.append(
                f'READS count too low ({min_reads:.0f}) — adapter trimming '
                f'or alignment may have failed'
            )

    # Quantitative sanity: no single anticodon should dominate
    rpm_rows = [row for row in rows if row[0] not in ('READS', 'UndetNNN')]
    if rpm_rows:
        for col_idx in range(1, len(header)):
            col_vals = [float(r[col_idx]) for r in rpm_rows]
            total = sum(col_vals)
            if total > 0:
                max_frac = max(col_vals) / total
                if max_frac > 0.80:
                    max_ac = rpm_rows[col_vals.index(max(col_vals))][0]
                    failures.append(
                        f'{header[col_idx]}: {max_ac} has {max_frac:.0%} of total '
                        f'RPM — likely too few reads quantified'
                    )

        nonzero = sum(1 for r in rpm_rows if any(float(v) > 0 for v in r[1:]))
        if nonzero < 3:
            failures.append(
                f'Only {nonzero} anticodons detected (expected >= 3 even on mini genome)'
            )

    # K562-specific: check replicate correlation
    # With only 50K reads on a mini genome, we expect some anticodons to show up
    # but replicates should still be broadly similar
    if len(header) >= 3:
        k562_cols = {}
        for i, col in enumerate(header[1:], 1):
            for sample in K562_SAMPLES:
                if sample in col:
                    k562_cols[sample] = i

        if len(k562_cols) == 2:
            vals1, vals2 = [], []
            for row in rows:
                if row[0] == 'READS':
                    continue
                try:
                    v1 = float(row[list(k562_cols.values())[0]])
                    v2 = float(row[list(k562_cols.values())[1]])
                    if v1 > 0 or v2 > 0:
                        vals1.append(v1)
                        vals2.append(v2)
                except (ValueError, IndexError):
                    continue

            if len(vals1) >= 3:
                # Simple rank correlation check
                # (can't import scipy in a basic test, use manual Spearman)
                n = len(vals1)
                rank1 = [sorted(vals1).index(v) for v in vals1]
                rank2 = [sorted(vals2).index(v) for v in vals2]
                d_sq = sum((r1 - r2) ** 2 for r1, r2 in zip(rank1, rank2))
                rho = 1 - (6 * d_sq) / (n * (n**2 - 1))
                if rho < 0.5:
                    failures.append(
                        f'K562 replicate Spearman rho = {rho:.3f} (expected >= 0.5 '
                        f'even with 50K reads on mini genome)'
                    )

    # --- results_rpm_nomod.csv ---
    nomod_path = results_dir / 'results_rpm_nomod.csv'
    if not nomod_path.exists():
        failures.append('results_rpm_nomod.csv does not exist')
    elif nomod_path.stat().st_size == 0:
        failures.append('results_rpm_nomod.csv is empty')

    # --- multiqc_report.html ---
    multiqc_path = results_dir / 'multiqc_report.html'
    if not multiqc_path.exists():
        failures.append('multiqc_report.html does not exist')
    elif multiqc_path.stat().st_size < 1000:
        failures.append(
            f'multiqc_report.html too small ({multiqc_path.stat().st_size} bytes)'
        )

    # --- pipeline_info/trace.txt ---
    trace_path = results_dir / 'pipeline_info' / 'trace.txt'
    if trace_path.exists():
        with open(trace_path) as f:
            trace_text = f.read()
        for proc in [
            'TRNASCAN_SE', 'MASK_GENOME', 'BUILD_PRETRNA',
            'BUILD_ARTIFICIAL_GENOME', 'BUILD_MATURE_LIBRARY',
            'INDEX_REFERENCES', 'TRIM_READS', 'MAP_ARTIFICIAL',
            'MAP_CLUSTERS', 'FILTER_TOPSCORE', 'SALMON_QUANT',
            'SALMON_TO_RPM', 'CALL_VARIANTS', 'COUNT_ALLELES',
            'ADJUST_MODIFICATIONS', 'AGGREGATE_RESULTS', 'MULTIQC',
        ]:
            if proc not in trace_text:
                failures.append(f'Process {proc} not found in trace.txt')

        if 'FAILED' in trace_text:
            n_failed = trace_text.count('FAILED')
            failures.append(f'{n_failed} FAILED process(es) in trace.txt')
    else:
        failures.append('trace.txt does not exist')

    return failures


def main():
    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} <results_dir>', file=sys.stderr)
        sys.exit(2)

    results_dir = Path(sys.argv[1])
    if not results_dir.is_dir():
        print(f'ERROR: {results_dir} is not a directory', file=sys.stderr)
        sys.exit(2)

    failures = validate(results_dir)

    if failures:
        print(f'FAILED ({len(failures)} issue(s)):')
        for f in failures:
            print(f'  - {f}')
        sys.exit(1)
    else:
        print('All K562 Python validations passed.')
        sys.exit(0)


if __name__ == '__main__':
    main()
