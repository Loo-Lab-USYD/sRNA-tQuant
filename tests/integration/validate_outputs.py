#!/usr/bin/env python3
"""
Validate integration test outputs from the tRNA mapping pipeline.

Usage:
    python3 validate_outputs.py <results_dir>

Checks RPM matrix structure, anticodon names, value ranges, and output
file integrity. Returns exit code 0 if all checks pass, 1 otherwise.
"""

import csv
import re
import sys
from pathlib import Path


# Known anticodon pattern: optional 'i' prefix, amino acid name, 3-letter codon
# e.g. AlaCGC, iMetCAT, SeC(e)TCA, SupTCA
ANTICODON_PATTERN = re.compile(r'^i?[A-Z][a-zA-Z()]+[ACTGN]{3}$')


def validate(results_dir: Path) -> list[str]:
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

    # Check sample IDs in header (may have _adjusted suffix from mod calling)
    sample_ids = header[1:]  # first column is index
    for expected in ['HCT116_test', 'HeLa_test']:
        if not any(expected in s for s in sample_ids):
            failures.append(f'Sample {expected} not in RPM header: {sample_ids}')

    # Check anticodon names in index
    index_values = [row[0] for row in rows]
    anticodon_rows = [v for v in index_values if v != 'READS']

    if 'READS' not in index_values:
        failures.append('READS row not found in RPM matrix')

    bad_names = []
    for name in anticodon_rows:
        if not ANTICODON_PATTERN.match(name):
            bad_names.append(name)
    if bad_names:
        failures.append(
            f'Invalid anticodon names (possible bug #12 regression): {bad_names[:5]}'
        )

    # Check no raw FASTA headers leaked (::chr pattern)
    for name in index_values:
        if '::' in name or '>' in name:
            failures.append(f'Raw FASTA header in index: {name}')
            break

    # Check RPM values are non-negative numbers
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
        # Sanity: even on mini genome with 50K reads, expect > 10 tRNA reads
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

        # At least a few anticodons should be detected
        nonzero = sum(1 for r in rpm_rows if any(float(v) > 0 for v in r[1:]))
        if nonzero < 3:
            failures.append(
                f'Only {nonzero} anticodons detected (expected >= 3 even on mini genome)'
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

    # --- pipeline_info ---
    trace_path = results_dir / 'pipeline_info' / 'trace.txt'
    if trace_path.exists():
        with open(trace_path) as f:
            trace_lines = f.readlines()
        # Check all key processes appear in trace
        trace_text = ''.join(trace_lines)
        required_procs = [
            'TRNASCAN_SE', 'MASK_GENOME', 'BUILD_PRETRNA',
            'BUILD_ARTIFICIAL_GENOME', 'BUILD_MATURE_LIBRARY',
            'INDEX_REFERENCES', 'TRIM_READS', 'MAP_ARTIFICIAL',
            'MAP_CLUSTERS', 'FILTER_TOPSCORE', 'SALMON_QUANT',
            'SALMON_TO_RPM', 'AGGREGATE_RESULTS', 'MULTIQC',
        ]
        # Mod calling processes are optional (skip_modification_calling)
        optional_procs = ['CALL_VARIANTS', 'COUNT_ALLELES', 'ADJUST_MODIFICATIONS']
        for proc in required_procs:
            if proc not in trace_text:
                failures.append(f'Process {proc} not found in trace.txt')
        # If any optional mod process ran, all must have ran
        mod_found = [p for p in optional_procs if p in trace_text]
        if 0 < len(mod_found) < len(optional_procs):
            missing = [p for p in optional_procs if p not in trace_text]
            failures.append(f'Partial mod calling: missing {missing}')

        # Check no FAILED status in trace
        failed = [l for l in trace_lines if 'FAILED' in l]
        if failed:
            failures.append(f'{len(failed)} FAILED process(es) in trace.txt')
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
        print('All Python validations passed.')
        sys.exit(0)


if __name__ == '__main__':
    main()
