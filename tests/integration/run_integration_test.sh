#!/usr/bin/env bash
#
# Integration test: run the full tRNA mapping pipeline on a subsampled dataset
#
# Uses chrM + chr6_GL000250v2_alt (4.8 MB, ~43 tRNAs) instead of full hg38 (3.3 GB)
# and 50K reads from 2 samples. Should complete in ~1-3 minutes.
#
# Usage:
#   bash tests/integration/run_integration_test.sh [--clean] [--keep-work]
#
# Options:
#   --clean       Remove work dir and results before running
#   --keep-work   Don't remove work dir after successful run (for debugging)

set -uo pipefail

# Java 21 from mamba env (system Java is too new for Nextflow)
export JAVA_HOME="/home/labhund/mamba_envs/trna_validation/lib/jvm"
export PATH="/home/labhund/mamba_envs/trna_validation/bin:$PATH"
export TERM="${TERM:-xterm}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
FIXTURES_DIR="${SCRIPT_DIR}/fixtures"
RESULTS_DIR="${SCRIPT_DIR}/results"
WORK_DIR="${SCRIPT_DIR}/work"
LOG_FILE="${SCRIPT_DIR}/pipeline.log"

CLEAN=false
KEEP_WORK=false
for arg in "$@"; do
    case "$arg" in
        --clean) CLEAN=true ;;
        --keep-work) KEEP_WORK=true ;;
    esac
done

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

pass() { echo -e "  ${GREEN}PASS${NC} $1"; }
fail() { echo -e "  ${RED}FAIL${NC} $1"; FAILURES=$((FAILURES + 1)); }
info() { echo -e "${YELLOW}>>>${NC} $1"; }

FAILURES=0

# ── Pre-flight checks ──────────────────────────────────────────────────

info "Checking fixtures..."

if [[ ! -f "${FIXTURES_DIR}/mini_genome.fa" ]]; then
    echo "ERROR: ${FIXTURES_DIR}/mini_genome.fa not found."
    echo "Run setup_fixtures.sh first, or see README."
    exit 1
fi

for f in HCT116_test.fastq.gz HeLa_test.fastq.gz; do
    if [[ ! -f "${FIXTURES_DIR}/${f}" ]]; then
        echo "ERROR: ${FIXTURES_DIR}/${f} not found."
        exit 1
    fi
done

# Resolve absolute paths in samplesheet
SAMPLESHEET="${FIXTURES_DIR}/samplesheet_resolved.csv"
sed "s|FIXTURES_DIR|${FIXTURES_DIR}|g" "${FIXTURES_DIR}/samplesheet.csv" > "${SAMPLESHEET}"

if $CLEAN; then
    info "Cleaning previous run..."
    rm -rf "${WORK_DIR}" "${RESULTS_DIR}" "${LOG_FILE}"
fi

# ── Run pipeline ────────────────────────────────────────────────────────

info "Running pipeline (mini genome: chrM + chr6_GL000250v2_alt, 2 samples x 50K reads)..."

START_TIME=$(date +%s)

set +e
nextflow run "${PROJECT_DIR}/pipeline/main.nf" \
    --input "${SAMPLESHEET}" \
    --genome_fasta "${FIXTURES_DIR}/mini_genome.fa" \
    --outdir "${RESULTS_DIR}" \
    --max_cpus 4 \
    --max_memory '8.GB' \
    --bowtie2_seed_length 14 \
    --bowtie2_seed_extensions 25 \
    -profile docker \
    -work-dir "${WORK_DIR}" \
    2>&1 | tee "${LOG_FILE}"

NF_EXIT=${PIPESTATUS[0]}
set -e
END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))

echo ""
info "Pipeline finished in ${ELAPSED}s (exit code ${NF_EXIT})"
echo ""

# ── Assertions ──────────────────────────────────────────────────────────

info "Running assertions..."

# 1. Nextflow exit code
if [[ $NF_EXIT -eq 0 ]]; then
    pass "Nextflow exit code 0"
else
    fail "Nextflow exit code ${NF_EXIT}"
fi

# 2. All expected processes completed (check the log)
EXPECTED_PROCESSES=(
    "TRNASCAN_SE"
    "MASK_GENOME"
    "BUILD_PRETRNA"
    "BUILD_ARTIFICIAL_GENOME"
    "BUILD_MATURE_LIBRARY"
    "INDEX_REFERENCES"
    "TRIM_READS"
    "MAP_ARTIFICIAL"
    "MAP_CLUSTERS"
    "FILTER_TOPSCORE"
    "SALMON_QUANT"
    "SALMON_TO_RPM"
    "AGGREGATE_RESULTS"
    "MULTIQC"
)

# Modification calling processes (only expected when not skipped)
MOD_PROCESSES=(
    "CALL_VARIANTS"
    "COUNT_ALLELES"
    "ADJUST_MODIFICATIONS"
)

for proc in "${EXPECTED_PROCESSES[@]}"; do
    if grep -q "${proc}" "${LOG_FILE}"; then
        pass "Process ${proc} ran"
    else
        fail "Process ${proc} not found in log"
    fi
done

# Check mod calling processes: if any ran, all must have ran
MOD_FOUND=0
for proc in "${MOD_PROCESSES[@]}"; do
    if grep -q "${proc}" "${LOG_FILE}"; then
        MOD_FOUND=$((MOD_FOUND + 1))
        pass "Process ${proc} ran (mod calling enabled)"
    fi
done
if [[ $MOD_FOUND -gt 0 && $MOD_FOUND -lt ${#MOD_PROCESSES[@]} ]]; then
    fail "Partial mod calling: only ${MOD_FOUND}/${#MOD_PROCESSES[@]} mod processes ran"
elif [[ $MOD_FOUND -eq 0 ]]; then
    pass "Mod calling skipped (skip_modification_calling=true)"
fi

# 3. No FAILED processes in log
FAILED_COUNT=$(grep -c 'FAILED\|Error executing process' "${LOG_FILE}" || true)
if [[ $FAILED_COUNT -eq 0 ]]; then
    pass "No FAILED processes"
else
    fail "${FAILED_COUNT} FAILED process(es) found in log"
fi

# 4. Output files exist and are non-empty
for f in results_rpm.csv results_rpm_nomod.csv; do
    if [[ -s "${RESULTS_DIR}/${f}" ]]; then
        pass "${f} exists and is non-empty"
    else
        fail "${f} missing or empty"
    fi
done

if [[ -f "${RESULTS_DIR}/multiqc_report.html" ]]; then
    pass "multiqc_report.html exists"
else
    fail "multiqc_report.html missing"
fi

# 5. RPM matrix content checks
if [[ -s "${RESULTS_DIR}/results_rpm.csv" ]]; then
    # Check column count (should have sample columns)
    NCOLS=$(head -1 "${RESULTS_DIR}/results_rpm.csv" | awk -F',' '{print NF}')
    if [[ $NCOLS -ge 3 ]]; then
        pass "RPM matrix has ${NCOLS} columns (>= 3: index + 2 samples)"
    else
        fail "RPM matrix has only ${NCOLS} columns (expected >= 3)"
    fi

    # Check row count (should have anticodon rows)
    NROWS=$(wc -l < "${RESULTS_DIR}/results_rpm.csv")
    if [[ $NROWS -ge 5 ]]; then
        pass "RPM matrix has ${NROWS} rows (>= 5)"
    else
        fail "RPM matrix has only ${NROWS} rows (expected >= 5)"
    fi

    # Check that headers contain sample IDs
    HEADER=$(head -1 "${RESULTS_DIR}/results_rpm.csv")
    if echo "$HEADER" | grep -q "HCT116_test" && echo "$HEADER" | grep -q "HeLa_test"; then
        pass "RPM matrix headers contain both sample IDs"
    else
        fail "RPM matrix headers missing sample IDs: ${HEADER}"
    fi

    # Check no raw FASTA headers leaked into index (bug #12 regression)
    if grep -q '::chr' "${RESULTS_DIR}/results_rpm.csv"; then
        fail "RPM matrix contains raw FASTA headers (::chr) — bug #12 regression"
    else
        pass "RPM matrix index has clean anticodon names (no ::chr)"
    fi
fi

# 6. Quantitative sanity: check fastp actually trimmed adapters
for sample in HCT116_test HeLa_test; do
    FASTP_JSON=$(find "${WORK_DIR}" -name "${sample}_fastp.json" -type f 2>/dev/null | head -1)
    if [[ -n "$FASTP_JSON" && -f "$FASTP_JSON" ]]; then
        TRIM_RESULT=$(python3 -c "
import json
with open('$FASTP_JSON') as f:
    d = json.load(f)
ac = d.get('adapter_cutting', {})
trimmed = ac.get('adapter_trimmed_reads', 0)
total = d['summary']['before_filtering']['total_reads']
pct = 100 * trimmed / total if total > 0 else 0
print(f'{pct:.1f} {trimmed} {total}')
" 2>/dev/null || echo "0.0 0 0")
        TRIM_PCT=$(echo "$TRIM_RESULT" | awk '{print $1}')
        TRIMMED_READS=$(echo "$TRIM_RESULT" | awk '{print $2}')
        TOTAL_READS=$(echo "$TRIM_RESULT" | awk '{print $3}')
        # At least 1% of reads should have adapter (even on mini genome subsets)
        if python3 -c "exit(0 if $TRIM_PCT > 1.0 else 1)"; then
            pass "fastp ${sample}: ${TRIM_PCT}% adapter-trimmed (${TRIMMED_READS}/${TOTAL_READS})"
        else
            fail "fastp ${sample}: only ${TRIM_PCT}% adapter-trimmed — wrong adapter or untrimmed data"
        fi
    fi
done

# 7. Run detailed Python validation if pytest is available
if command -v python3 &>/dev/null; then
    info "Running Python validation..."
    if python3 "${SCRIPT_DIR}/validate_outputs.py" "${RESULTS_DIR}"; then
        pass "Python validation passed"
    else
        fail "Python validation failed"
    fi
fi

# ── Summary ─────────────────────────────────────────────────────────────

echo ""
if [[ $FAILURES -eq 0 ]]; then
    echo -e "${GREEN}All assertions passed!${NC} (${ELAPSED}s)"
else
    echo -e "${RED}${FAILURES} assertion(s) failed.${NC} (${ELAPSED}s)"
fi

# Cleanup
if ! $KEEP_WORK && [[ $FAILURES -eq 0 ]]; then
    info "Cleaning work directory..."
    rm -rf "${WORK_DIR}"
fi

exit $FAILURES
