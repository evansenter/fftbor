# FFTbor Modernization TODO

Remaining modernization opportunities for the FFTbor codebase.

## High Priority

### 1. Fix unused parameter warnings
**File:** `delta.cpp:402-405`

The `populate_matrices()` function has 6 unused parameters:
- `Z`, `ZB`, `ZM`, `ZM1`, `solutions`, `sequence_length`

Either remove these parameters or document why they're needed for future use.

### 2. Remove dead legacy arrays
**File:** `globals.cpp`

The following arrays (~50 lines) were only needed for ViennaRNA library compatibility and are now unused:
- `stack37`, `enthalpies`, `entropies`
- `hairpin37`, `bulge37`, `internal_loop37`
- `old_mismatch_37`, `mismatchI37`, `mismatchH37`, `mismatchM37`, `mism_H`
- `dangle5_37`, `dangle3_37`, `dangle3_H`, `dangle5_H`
- `int11_37`, `int11_H`, `int21_37`, `int21_H`, `int22_37`, `int22_H`
- `ML_BASE37`, `ML_closing37`, `ML_intern37`
- `F_ninio37`, `TerminalAU`, `DuplexInit`
- `Tetraloops`, `TETRA_ENERGY37`, `TETRA_ENTH37`
- `Triloops`, `Triloop_E37`

### 3. Full Context migration
**Files:** `delta.cpp`, `main.cpp`, `globals.cpp`

Replace global variables with `fftbor::Context` passed through functions:
- Migrate `N`, `PRECISION`, `WINDOW_SIZE`, `MIN_WINDOW_SIZE`, `temperature`, `ENERGY`, `P`
- Update function signatures to accept `const Context&` or `Context&`
- Delete `globals.cpp` after migration complete

## Medium Priority

### 4. Modern string APIs
**Files:** `delta.h`, `misc.h`

Replace C-style string parameters with `std::string_view`:
- `neighbours(const char* input_sequence, ...)` → `neighbours(std::string_view input_sequence, ...)`
- `get_base_pair_list(const char* sec_str)` → `get_base_pair_list(std::string_view sec_str)`

### 5. RAII for FILE handles
**Files:** `parameter_parser.cpp`, `main.cpp`

Wrap `fopen`/`fclose` with RAII:
```cpp
struct FileCloser {
    void operator()(FILE* f) const { if (f) fclose(f); }
};
using FilePtr = std::unique_ptr<FILE, FileCloser>;
```

### 6. Replace printf with iostream
**Files:** `delta.cpp`, `main.cpp`

Replace C-style I/O with C++ streams for consistency:
- `printf(...)` → `std::cout << ...`
- `fprintf(stderr, ...)` → `std::cerr << ...`

## Lower Priority

### 7. C++ style casts
**Files:** Various

Replace C-style casts with explicit C++ casts:
- `(int*)xcalloc(...)` → `static_cast<int*>(xcalloc(...))`
- `(char*)xcalloc(...)` → `static_cast<char*>(xcalloc(...))`

### 8. Safer string functions
**Files:** `delta.cpp`, `main.cpp`

Replace potentially unsafe string functions:
- `strcpy` → `std::copy` or bounds-checked alternative
- `strncpy` → `std::string` or `std::copy_n`
- `sscanf` → `std::istringstream` or `std::from_chars`

---

## Completed Modernizations

- [x] Upgrade to C++20 (concepts support)
- [x] Replace `exit(1)` with `std::runtime_error` exceptions
- [x] Add `fftbor::Context` struct for state encapsulation
- [x] Add C++20 concepts: `SequenceView`, `BasePairSpan`
- [x] Convert macros to inline functions
- [x] Add named constants for magic numbers
- [x] Update `scale_parameters()` to take explicit temperature parameter
- [x] Add `load_parameters(Context&)` convenience function
- [x] Smart pointers for memory management
- [x] Bazel build system
- [x] GoogleTest test suite
