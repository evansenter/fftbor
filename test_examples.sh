#!/bin/bash
# Test script for FFTbor examples
# Runs various RNA structures through FFTbor and validates output
#
# Usage:
#   ./test_examples.sh        # Run tests in quiet mode
#   ./test_examples.sh -v     # Run tests in verbose mode (show full output)

set -e

VERBOSE=0
PARAM_FILE="rna_turner2004.par"
PASS=0
FAIL=0

# Parse command line arguments
while getopts "v" opt; do
    case $opt in
        v) VERBOSE=1 ;;
        *) echo "Usage: $0 [-v]"; exit 1 ;;
    esac
done

# Check if rebuild is needed
check_and_rebuild() {
    local dominated_sources=("*.cpp" "*.h" "BUILD.bazel" "MODULE.bazel")
    local needs_rebuild=0

    # Check if binary exists
    if [[ ! -f "./bazel-bin/FFTbor" ]]; then
        echo "Binary not found, building..."
        needs_rebuild=1
    else
        # Check if any source files are newer than the binary
        local binary_time=$(stat -f %m "./bazel-bin/FFTbor" 2>/dev/null || stat -c %Y "./bazel-bin/FFTbor" 2>/dev/null)

        for pattern in "${dominated_sources[@]}"; do
            for file in $pattern; do
                if [[ -f "$file" ]]; then
                    local file_time=$(stat -f %m "$file" 2>/dev/null || stat -c %Y "$file" 2>/dev/null)
                    if [[ "$file_time" -gt "$binary_time" ]]; then
                        echo "Source file $file has changed, rebuilding..."
                        needs_rebuild=1
                        break 2
                    fi
                fi
            done
        done
    fi

    if [[ $needs_rebuild -eq 1 ]]; then
        echo "Running: bazel build //:FFTbor"
        if ! bazel build //:FFTbor 2>&1; then
            echo "Error: Build failed"
            exit 1
        fi
        echo "Build successful"
        echo ""
    fi
}

# Rebuild if needed
check_and_rebuild

# Determine the FFTbor binary location
if [[ -f "./bazel-bin/FFTbor" ]]; then
    FFTBOR="./bazel-bin/FFTbor"
elif [[ -f "./bin/FFTbor" ]]; then
    FFTBOR="./bin/FFTbor"
else
    echo "Error: FFTbor binary not found."
    exit 1
fi

# Helper function to run a test (quiet mode)
run_test_quiet() {
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

# Helper function to run a test and show output (verbose mode)
run_test_verbose() {
    local name="$1"
    local seq="$2"
    local str="$3"
    local extra_args="${4:-}"

    echo "=== $name ==="
    if $FFTBOR -E "$PARAM_FILE" $extra_args "$seq" "$str"; then
        ((PASS++))
    else
        ((FAIL++))
    fi
    echo ""
}

# Wrapper that calls the appropriate function based on mode
run_test() {
    if [[ $VERBOSE -eq 1 ]]; then
        run_test_verbose "$@"
    else
        run_test_quiet "$@"
    fi
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
