#!/bin/bash
# Test script for FFTbor examples
# Runs various RNA structures through FFTbor and validates output

set -e

# Determine the FFTbor binary location
if [[ -f "./bazel-bin/FFTbor" ]]; then
    FFTBOR="./bazel-bin/FFTbor"
elif [[ -f "./bin/FFTbor" ]]; then
    FFTBOR="./bin/FFTbor"
else
    echo "Error: FFTbor binary not found. Build with 'bazel build //:FFTbor' first."
    exit 1
fi

PARAM_FILE="rna_turner2004.par"
PASS=0
FAIL=0

# Helper function to run a test
run_test() {
    local name="$1"
    local seq="$2"
    local str="$3"
    local extra_args="${4:-}"

    echo -n "Testing $name... "

    if output=$($FFTBOR -E "$PARAM_FILE" $extra_args "$seq" "$str" 2>&1); then
        # Check that output contains expected format (k and p(k) columns)
        if echo "$output" | grep -q "^k	p(k)$"; then
            echo "PASS"
            ((PASS++))
            return 0
        else
            echo "FAIL (unexpected output format)"
            ((FAIL++))
            return 1
        fi
    else
        echo "FAIL (exit code $?)"
        echo "$output"
        ((FAIL++))
        return 1
    fi
}

# Helper function to run a test and show output
run_test_verbose() {
    local name="$1"
    local seq="$2"
    local str="$3"
    local extra_args="${4:-}"

    echo "=== $name ==="
    $FFTBOR -E "$PARAM_FILE" $extra_args "$seq" "$str"
    echo ""
}

echo "FFTbor Example Test Suite"
echo "========================="
echo "Using binary: $FFTBOR"
echo ""

# Test 1: Simple hairpin
run_test "Simple hairpin (12nt)" \
    "GCGCAAAAGCGC" \
    "((((....))))"

# Test 2: Two-stem structure
run_test "Two-stem structure (28nt)" \
    "GGGGAAAACCCCUUUUGGGGAAAACCCC" \
    "((((....))))....((((....))))"

# Test 3: tRNA-like cloverleaf (73nt)
run_test "tRNA cloverleaf (73nt)" \
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA" \
    "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))."

# Test 4: Multiloop structure (23nt)
run_test "Multiloop structure (23nt)" \
    "GGAAUCCGAAAUUGGAAUUCCGG" \
    "((..(((....)))..((.))))"

# Test 5: Temperature variation (25°C)
run_test "tRNA at 25°C" \
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA" \
    "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))." \
    "-T 25"

# Test 6: Temperature variation (50°C)
run_test "tRNA at 50°C" \
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA" \
    "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))." \
    "-T 50"

# Test 7: High precision output
run_test "tRNA high precision (P=8)" \
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA" \
    "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))." \
    "-P 8"

echo ""
echo "========================="
echo "Results: $PASS passed, $FAIL failed"

if [[ $FAIL -eq 0 ]]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed."
    exit 1
fi
