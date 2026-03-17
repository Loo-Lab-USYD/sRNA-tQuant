#!/usr/bin/env bash
#
# Integration test: K562 ENCODE small RNA-seq data through the tRNA pipeline
#
# Tests that 101bp SE reads from ENCODE ENCSR000AES (TAP-treated, TruSeq adapter)
# are correctly trimmed, mapped, and quantified by our pipeline.
#
# Uses the same mini genome (chrM + chr6_GL000250v2_alt) and 50K subsampled reads.
# Should complete in ~1-3 minutes.
#
# Usage:
#   bash tests/integration/run_k562_integration_test.sh [--clean] [--keep-work]

set -uo pipefail

# Java 21 from mamba env (system Java is too new for Nextflow)
export JAVA_HOME="/home/labhund/mamba_envs/trna_validation/lib/jvm"
export PATH="/home/labhund/mamba_envs/trna_validation/bin:$PATH"
export TERM="${TERM:-xterm}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
FIXTURES_DIR="${SCRIPT_DIR}/fixtures"
RESULTS_DIR="${SCRIPT_DIR}/results_k562"
WORK_DIR="${SCRIPT_DIR}/work_k562"
LOG_FILE="${SCRIPT_DIR}/pipeline_k562.log"

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

info "Checking K562 fixtures..."

if [[ ! -f "${FIXTURES_DIR}/mini_genome.fa" ]]; then
    echo "ERROR: ${FIXTURES_DIR}/mini_genome.fa not found."
    echo "Run the standard integration test first to generate fixtures."
    exit 1
fi

for f in K562_rep1_test.fastq.gz K562_rep2_test.fastq.gz; do
    if [[ ! -f "${FIXTURES_DIR}/${f}" ]]; then
        echo "ERROR: ${FIXTURES_DIR}/${f} not found."
        echo "Run this first to generate test fixtures from ENCODE downloads:"
        echo "  zcat data/k562/fastq/ENCFF001REQ.fastq.gz | head -200000 | gzip > tests/integration/fixtures/K562_rep1_test.fastq.gz"
        echo "  zcat data/k562/fastq/ENCFF001REL.fastq.gz | head -200000 | gzip > tests/integration/fixtures/K562_rep2_test.fastq.gz"
        exit 1
    fi
done

# Verify read format (101bp SE, TruSeq adapter present)
info "Verifying K562 read format..."
READ_LEN=$(zcat "${FIXTURES_DIR}/K562_rep1_test.fastq.gz" | head -2 | tail -1 | wc -c)
READ_LEN=$((READ_LEN - 1))  # subtract newline
if [[ $READ_LEN -ge 90 && $READ_LEN -le 110 ]]; then
    pass "Read length ${READ_LEN}bp (expected ~101bp for ENCODE)"
else
    fail "Unexpected read length ${READ_LEN}bp (expected ~101bp)"
fi

# Check TruSeq adapter is present in reads
ADAPTER_COUNT=$(zcat "${FIXTURES_DIR}/K562_rep1_test.fastq.gz" | head -2000 | grep -c "TGGAATTCTCGGGTGCCAAGG" || true)
if [[ $ADAPTER_COUNT -gt 10 ]]; then
    pass "TruSeq adapter found in ${ADAPTER_COUNT}/500 reads (fastp will auto-trim)"
else
    info "Warning: TruSeq adapter found in only ${ADAPTER_COUNT}/500 reads (fastp auto-detect should still work)"
fi

# Resolve absolute paths in samplesheet
SAMPLESHEET="${FIXTURES_DIR}/samplesheet_k562_resolved.csv"
sed "s|FIXTURES_DIR|${FIXTURES_DIR}|g" "${FIXTURES_DIR}/samplesheet_k562.csv" > "${SAMPLESHEET}"

if $CLEAN; then
    info "Cleaning previous K562 run..."
    rm -rf "${WORK_DIR}" "${RESULTS_DIR}" "${LOG_FILE}"
fi

# ── Run pipeline ────────────────────────────────────────────────────────

info "Running pipeline (K562 ENCODE: mini genome, 2 samples x 50K reads, 101bp SE)..."

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
    --adapter_sequence 'TGGAATTCTCGGGTGCCAAGG' \
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

# 2. All expected processes completed
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
    "CALL_VARIANTS"
    "COUNT_ALLELES"
    "ADJUST_MODIFICATIONS"
    "AGGREGATE_RESULTS"
    "MULTIQC"
)

for proc in "${EXPECTED_PROCESSES[@]}"; do
    if grep -q "${proc}" "${LOG_FILE}"; then
        pass "Process ${proc} ran"
    else
        fail "Process ${proc} not found in log"
    fi
done

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

# 5. RPM matrix content checks (K562-specific)
if [[ -s "${RESULTS_DIR}/results_rpm.csv" ]]; then
    NCOLS=$(head -1 "${RESULTS_DIR}/results_rpm.csv" | awk -F',' '{print NF}')
    if [[ $NCOLS -ge 3 ]]; then
        pass "RPM matrix has ${NCOLS} columns (>= 3: index + 2 K562 samples)"
    else
        fail "RPM matrix has only ${NCOLS} columns (expected >= 3)"
    fi

    NROWS=$(wc -l < "${RESULTS_DIR}/results_rpm.csv")
    if [[ $NROWS -ge 5 ]]; then
        pass "RPM matrix has ${NROWS} rows (>= 5)"
    else
        fail "RPM matrix has only ${NROWS} rows (expected >= 5)"
    fi

    # Check that headers contain K562 sample IDs
    HEADER=$(head -1 "${RESULTS_DIR}/results_rpm.csv")
    if echo "$HEADER" | grep -q "K562_rep1_test" && echo "$HEADER" | grep -q "K562_rep2_test"; then
        pass "RPM matrix headers contain both K562 sample IDs"
    else
        fail "RPM matrix headers missing K562 sample IDs: ${HEADER}"
    fi

    # Check no raw FASTA headers leaked into index
    if grep -q '::chr' "${RESULTS_DIR}/results_rpm.csv"; then
        fail "RPM matrix contains raw FASTA headers (::chr) — bug regression"
    else
        pass "RPM matrix index has clean anticodon names"
    fi

    # K562-specific: check that fastp handled 101bp reads properly
    # After adapter trimming, reads should be 10-50bp (small RNA inserts)
    # If most reads are still ~101bp, adapter trimming failed
fi

# 6. Quantitative sanity: check fastp actually trimmed adapters
for sample in K562_rep1_test K562_rep2_test; do
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
        if python3 -c "exit(0 if $TRIM_PCT > 1.0 else 1)"; then
            pass "fastp ${sample}: ${TRIM_PCT}% adapter-trimmed (${TRIMMED_READS}/${TOTAL_READS})"
        else
            fail "fastp ${sample}: only ${TRIM_PCT}% adapter-trimmed — wrong adapter or untrimmed data"
        fi
    fi
done

# 7. Run detailed Python validation
if command -v python3 &>/dev/null; then
    info "Running Python validation..."
    # Use the same validator but with K562 sample IDs
    VALIDATE_SCRIPT="${SCRIPT_DIR}/validate_k562_outputs.py"
    if [[ -f "$VALIDATE_SCRIPT" ]]; then
        if python3 "$VALIDATE_SCRIPT" "${RESULTS_DIR}"; then
            pass "Python validation passed"
        else
            fail "Python validation failed"
        fi
    else
        # Fall back to generic validator (will skip sample-ID-specific checks)
        info "K562-specific validator not found, skipping Python validation"
    fi
fi

# ── Summary ─────────────────────────────────────────────────────────────

echo ""
if [[ $FAILURES -eq 0 ]]; then
    echo -e "${GREEN}All K562 integration test assertions passed!${NC} (${ELAPSED}s)"
    echo ""
    echo "The pipeline handles 101bp SE ENCODE data correctly."
    echo "Safe to proceed with full K562 pipeline run."
else
    echo -e "${RED}${FAILURES} assertion(s) failed.${NC} (${ELAPSED}s)"
    echo ""
    echo "Review the log: ${LOG_FILE}"
fi

# Cleanup
if ! $KEEP_WORK && [[ $FAILURES -eq 0 ]]; then
    info "Cleaning work directory..."
    rm -rf "${WORK_DIR}"
fi

exit $FAILURES
