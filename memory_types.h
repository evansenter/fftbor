#ifndef MEMORY_TYPES_H
#define MEMORY_TYPES_H

#include <memory>
#include <vector>
#include <complex>
#include <cstdlib>
#include <string>
#include <string_view>
#include <span>
#include <concepts>
#include <cmath>

// Forward declaration
struct paramT;

namespace fftbor {

// Complex number type (matching delta.h)
using dcomplex = std::complex<double>;

// ============================================================================
// Custom Deleters for C-allocated memory
// ============================================================================

struct FreeDeleter {
    void operator()(void* ptr) const noexcept {
        std::free(ptr);
    }
};

// ============================================================================
// Smart Pointer Type Aliases
// ============================================================================

// Base pair list: allocated with xcalloc, needs free()
using BasePairListPtr = std::unique_ptr<int[], FreeDeleter>;

// Integer sequence: allocated with xcalloc, needs free()
using IntSequencePtr = std::unique_ptr<int[], FreeDeleter>;

// Character sequence: allocated with xcalloc, needs free()
using CharSequencePtr = std::unique_ptr<char[], FreeDeleter>;

// Wrap a raw char* (from xcalloc) in a smart pointer
inline CharSequencePtr make_char_sequence(char* raw_ptr) {
    return CharSequencePtr(raw_ptr);
}

// Parameter struct: allocated with calloc, needs free()
using ParamPtr = std::unique_ptr<paramT, FreeDeleter>;

// ============================================================================
// Matrix Types (replacing raw 2D/3D pointer arrays)
// ============================================================================

template<typename T>
using Matrix2D = std::vector<std::vector<T>>;

template<typename T>
using Matrix3D = std::vector<std::vector<std::vector<T>>>;

// Convenience aliases
using ComplexMatrix2D = Matrix2D<dcomplex>;
using ComplexMatrix3D = Matrix3D<dcomplex>;
using IntMatrix2D = Matrix2D<int>;

// ============================================================================
// Factory Functions
// ============================================================================

// Wrap a raw int* (from xcalloc) in a smart pointer
inline BasePairListPtr make_base_pair_list(int* raw_ptr) {
    return BasePairListPtr(raw_ptr);
}

// Wrap a raw paramT* (from calloc) in a smart pointer
inline ParamPtr make_param_ptr(paramT* raw_ptr) {
    return ParamPtr(raw_ptr);
}

// Create a 2D complex matrix initialized to zero
inline ComplexMatrix2D make_complex_matrix_2d(size_t rows, size_t cols) {
    return ComplexMatrix2D(rows, std::vector<dcomplex>(cols, dcomplex(0.0, 0.0)));
}

// Create a 2D int matrix initialized to a value
inline IntMatrix2D make_int_matrix_2d(size_t rows, size_t cols, int init_val = 0) {
    return IntMatrix2D(rows, std::vector<int>(cols, init_val));
}

// ============================================================================
// C++20 Concepts for Type Safety
// ============================================================================

template<typename T>
concept SequenceView = requires(T seq) {
    { seq.size() } -> std::convertible_to<size_t>;
    { seq[0] } -> std::convertible_to<char>;
    { seq.data() } -> std::convertible_to<const char*>;
};

template<typename T>
concept BasePairSpan = requires(T bpl) {
    { bpl.size() } -> std::convertible_to<size_t>;
    { bpl[0] } -> std::convertible_to<int>;
    { bpl.data() } -> std::convertible_to<const int*>;
};

// ============================================================================
// Physical Constants
// ============================================================================

// Gas constant R in kcal/(mol·K), scaled for internal energy units (1/10 cal/mol)
// R = 1.987204258 cal/(mol·K) = 0.001987204258 kcal/(mol·K)
// Scaled by 100 for internal units: R_scaled = 0.0019872... * 100 = 0.19872...
// Pre-multiplied with kelvin offset for efficiency
constexpr double GAS_CONSTANT_KCAL_PER_MOL_K = 0.0019872370936902486;

// Kelvin offset for converting Celsius to Kelvin
constexpr double KELVIN_OFFSET = 273.15;

// ============================================================================
// Parameter File Constants
// ============================================================================

constexpr int MAX_TETRALOOPS = 200;
constexpr int TETRALOOP_SEQ_LENGTH = 6;  // nucleotides per tetraloop sequence
constexpr int TETRALOOP_SEQ_SIZE = MAX_TETRALOOPS * (TETRALOOP_SEQ_LENGTH + 1) + 1;  // +1 for null terminator per seq, +1 for array

constexpr int MAX_TRILOOPS = 40;
constexpr int TRILOOP_SEQ_LENGTH = 5;    // nucleotides per triloop sequence
constexpr int TRILOOP_SEQ_SIZE = MAX_TRILOOPS * (TRILOOP_SEQ_LENGTH + 1) + 1;

constexpr int MAX_HEXALOOPS = 40;
constexpr int HEXALOOP_SEQ_LENGTH = 8;   // nucleotides per hexaloop sequence
constexpr int HEXALOOP_SEQ_SIZE = MAX_HEXALOOPS * (HEXALOOP_SEQ_LENGTH + 1) + 1;

// ============================================================================
// Algorithm Constants
// ============================================================================

constexpr int MIN_PAIR_DIST = 3;
constexpr int MAX_INTERIOR_DIST = 30;

// ============================================================================
// Context - Encapsulates all algorithm state (replaces globals)
// ============================================================================

struct Context {
    // Configuration (set by user/CLI)
    int precision = 4;
    int window_size = 0;
    int min_window_size = 0;
    double temperature = 37.0;
    std::string energy_file = "rna_turner2004.par";

    // Computed state
    int sequence_length = 0;
    double rt = 0.0;  // Pre-computed RT = GAS_CONSTANT_KCAL_PER_MOL_K * (T + KELVIN_OFFSET) * 100

    // Energy parameters (owned)
    ParamPtr params = nullptr;

    // Compute RT from temperature (in Celsius)
    // RT = R * T_kelvin, scaled by 100 for internal energy units
    void compute_rt() {
        rt = GAS_CONSTANT_KCAL_PER_MOL_K * (temperature + KELVIN_OFFSET) * 100;
    }

    // Get number of windows
    int num_windows() const {
        return window_size - min_window_size + 1;
    }

    // Get window size at index i
    int window_size_at(int i) const {
        return min_window_size + i;
    }
};

// ============================================================================
// Inline Helper Functions (replacing macros)
// ============================================================================

// Get root of unity raised to power: roots[(i * pow) % (n + 1)]
inline dcomplex root_pow(int i, int pow, int n, const std::vector<dcomplex>& roots) {
    return roots[(i * pow) % (n + 1)];
}

} // namespace fftbor

#endif // MEMORY_TYPES_H
