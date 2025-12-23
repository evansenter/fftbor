#ifndef MEMORY_TYPES_H
#define MEMORY_TYPES_H

#include <memory>
#include <vector>
#include <complex>
#include <cstdlib>

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

} // namespace fftbor

#endif // MEMORY_TYPES_H
