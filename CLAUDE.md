# CLAUDE.md - AI Assistant Guide for FFTbor

## Project Overview

FFTbor is a computational biology tool that computes the number of delta-neighbours and the partition function of delta-neighbours for RNA secondary structures. It uses Fast Fourier Transform (FFT) algorithms combined with thermodynamic energy calculations based on the Turner 2004 energy parameters.

**License:** Creative Commons Non-Commercial Share-Alike

## Repository Structure

```
fftbor/
├── main.cpp              # Entry point, CLI argument parsing
├── delta.cpp             # Core algorithm implementation
├── delta.h               # Algorithm interface/declarations
├── misc.cpp              # Utility functions (memory, parsing)
├── misc.h                # Utility function declarations
├── parameter_parser.cpp  # Native ViennaRNA parameter file parser
├── parameter_parser.h    # Parameter parser declarations
├── globals.cpp           # Global variable definitions
├── params.h              # Parameter structures
├── energy_const.h        # Energy calculation constants
├── energy_par.h          # Energy parameter declarations
├── rna_turner1999.par    # Turner 1999 energy parameters
├── rna_turner2004.par    # Turner 2004 energy parameters (default)
├── CMakeLists.txt        # CMake build configuration
├── cmake_modules/        # Custom CMake modules
│   └── FindFFTW.cmake    # FFTW library finder
├── tests/                # GoogleTest test suite
│   ├── test_main.cpp
│   ├── test_misc.cpp
│   ├── test_parameter_parser.cpp
│   ├── test_energy.cpp
│   └── test_integration.cpp
├── .github/workflows/    # CI/CD configuration
│   └── ci.yml
├── BUILD.bazel           # Bazel build configuration
├── MODULE.bazel          # Bazel module definition
├── TODO.md               # Remaining modernization tasks
└── build/                # Build output directory (CMake)
```

## Building the Project

### Dependencies

1. **C++ compiler** with OpenMP support (C++20 standard)
2. **CMake** (>= 3.16)
3. **FFTW** (= 3.3.x) - Fast Fourier Transform library

### Build with Bazel (Recommended)

```bash
bazel build //...
bazel test //tests:fftbor_tests
./bazel-bin/FFTbor -E rna_turner2004.par "GCGCAAAAGCGC" "((((....))))"
```

### Build with CMake

```bash
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
make install
```

### Build Configurations (CMake)

- **Default:** Optimized build with `-Ofast`, `-funroll-loops`, `-ffast-math`
- **Debug:** `cmake -DCMAKE_BUILD_TYPE=Debug ..` (enables `-g`)
- **Safe:** `cmake -DCMAKE_BUILD_TYPE=Safe ..` (uses `-O2` only)

## Key Source Files

### main.cpp (Entry Point)
- `main()`: Program entry, calls `read_input()` and `neighbours()`
- `read_input()`: Parses CLI options (`-T` temperature, `-E` energy file, `-P` precision, `-W` window size, `-M` min window size)
- `usage()`: Prints help information

### delta.cpp (Core Algorithm)
- `neighbours()`: Main algorithm entry point
  - Initializes energy parameters
  - Manages roots of unity for FFT computation
  - Calls `evaluate_z()` for each root of unity
  - Performs inverse DFT via `solve_system()`
- `evaluate_z()`: Computes partition function matrices (Z, ZB, ZM, ZM1)
- `hairpin_loop()`: Calculates hairpin loop energy
- `interior_loop()`: Calculates interior loop/bulge/stack energy
- `solve_system()`: Uses FFTW for DFT computation

### misc.cpp (Utilities)
- `xcalloc()`: Safe memory allocation wrapper
- `get_base_pair_list()`: Converts secondary structure notation to base pair list

## Code Conventions

### Naming
- **Global variables:** UPPERCASE (e.g., `N`, `PRECISION`, `WINDOW_SIZE`)
- **Functions:** snake_case (e.g., `read_input()`, `get_base_pair_list()`)
- **Local variables:** snake_case (e.g., `bp_list`, `sequence_length`, `run_length`)
- **Constants:** UPPERCASE with underscore (e.g., `MAX_INTERIOR_DIST`)
- **Types/structs:** Suffix with `T` (e.g., `paramT`, `model_detailsT`)

### Memory Management
- Smart pointers used for automatic memory management (RAII pattern)
- `xcalloc()` wrapper still available for allocation
- Custom smart pointer types in `memory_types.h` (e.g., `CharSequencePtr`, `BasePairListPtr`)

### Error Handling
- Exception-based error handling using `std::runtime_error`
- All `exit(1)` calls replaced with exceptions
- `main()` catches and handles exceptions gracefully

### Modern C++ Features
- **Context struct:** `fftbor::Context` in `memory_types.h` for state encapsulation
- **Concepts:** `SequenceView`, `BasePairSpan` for type safety
- **Inline functions:** `window_size_at()`, `num_windows()`, `fftbor::root_pow()` replace macros
- **Named constants:** `MIN_PAIR_DIST`, `MAX_INTERIOR_DIST`, etc. in `fftbor` namespace

### Data Representations
- **Base pair list:** Array where `bp_list[i]` = position paired with i (-1 if unpaired)
- **Sequence encoding:** 1=A, 2=C, 3=G, 4=U, 0=undefined
- **Complex numbers:** `std::complex<double>` (typedef `dcomplex`)

### Key Constants
- `MIN_PAIR_DIST = 3` - Minimum unpaired bases between paired positions
- `MAX_INTERIOR_DIST = 30` - Maximum unpaired bases in interior loops
- `MAXLOOP = 30` - Maximum single loop size
- `NBPAIRS = 7` - Number of base pair types

## Algorithm Overview

The FFTbor algorithm uses FFT to efficiently compute partition function coefficients:

1. **Roots of Unity:** Computes partition function at multiple roots of unity (e^(2πik/(n+1)))
2. **Symmetry Optimization:** Only computes ceil(run_length/2) roots due to complex conjugate symmetry
3. **Dynamic Programming:** Maintains four matrices:
   - `Z[i][j]`: Partition function for unpaired region [i,j]
   - `ZB[i][j]`: Partition function for region [i,j] with (i,j) base-paired
   - `ZM[i][j]`: Partition function for multiloop region [i,j]
   - `ZM1[i][j]`: Partition function for single component in multiloop [i,j]
4. **Inverse DFT:** Uses FFTW to transform point-value solutions back to coefficients

### Energy Calculations
- **Temperature scaling:** RT = 0.01 * 1.987 * (T + 273.15)
- **Boltzmann weight:** exp(-energy / RT)
- **Valid base pairs:** A-U, U-A, C-G, G-C, G-U, U-G

## CLI Usage

```bash
FFTbor [options] <sequence> <structure>
# or
FFTbor [options] < input.fasta
```

### Options
- `-T <temp>`: Temperature in Celsius (default: 37)
- `-E <file>`: Energy parameter file path
- `-P <precision>`: Decimal precision 0-9 (default: 4)
- `-W <size>`: Window size
- `-M <size>`: Minimum window size

## Debugging

Compile-time flags:
- `FFTBOR_DEBUG`: Enables verbose debug output with progress dots
- `ENERGY_DEBUG`: Enables energy calculation tracing
- `STRUCTURE_COUNT`: Enables structure counting feature

## Guidelines for AI Assistants

### When Modifying Code
1. Use modern C++20 features where appropriate (concepts, ranges, std::span, etc.)
2. Use the established naming conventions
3. Prefer RAII and smart pointers for new code
4. Keep energy calculation precision in mind
5. Test with both Turner 1999 and 2004 parameters

### When Adding Features
1. Add new functions to appropriate files (delta.cpp for algorithm, misc.cpp for utilities)
2. Declare new functions in corresponding header files
3. Update CMakeLists.txt if adding new source files
4. Ensure OpenMP compatibility for parallelization

### Common Tasks
- **Build:** `cd build && cmake .. && make`
- **Run tests:** `cd build && ctest` or `./bin/fftbor_tests`
- **Clean build:** `cd build && rm -rf *`

### Energy Parameters
- Default: `rna_turner2004.par`
- Alternative: `rna_turner1999.par`
- Format: ViennaRNA 2.x parameter format
- Installed to `bin/` directory alongside executable

### Important Notes
- VLAs replaced with `std::vector` for standard C++ compliance
- OpenMP is optional but recommended for parallel execution
- FFTW library handles the core FFT computations
- Native parameter parser included (no external ViennaRNA dependency)
- GoogleTest is used for unit and integration testing
