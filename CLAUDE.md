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
├── params.h              # Parameter structures (from ViennaRNA)
├── energy_const.h        # Energy calculation constants
├── energy_par.h          # Energy parameter declarations
├── rna_turner1999.par    # Turner 1999 energy parameters
├── rna_turner2004.par    # Turner 2004 energy parameters (default)
├── CMakeLists.txt        # CMake build configuration
├── cmake_modules/        # Custom CMake modules
│   └── FindFFTW.cmake    # FFTW library finder
└── build/                # Build output directory
```

## Building the Project

### Dependencies

1. **C++ compiler** with OpenMP support (C++98 standard)
2. **CMake** (>= 2.6)
3. **FFTW** (= 3.3.x) - Fast Fourier Transform library
4. **libRNA.a** from ViennaRNA package (= 2.x)

### Build Commands

```bash
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
make install
```

If `libRNA.a` is in a non-standard location, use `-DCMAKE_LIBRARY_PATH=/path/to/lib`.

### Build Configurations

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
  - Calls `evaluateZ()` for each root of unity
  - Performs inverse DFT via `solveSystem()`
- `evaluateZ()`: Computes partition function matrices (Z, ZB, ZM, ZM1)
- `hairpinloop()`: Calculates hairpin loop energy
- `interiorloop()`: Calculates interior loop/bulge/stack energy
- `solveSystem()`: Uses FFTW for DFT computation

### misc.cpp (Utilities)
- `xcalloc()`: Safe memory allocation wrapper
- `getBasePairList()`: Converts secondary structure notation to base pair list
- `min2()`, `max2()`: Comparison utilities

## Code Conventions

### Naming
- **Global variables:** UPPERCASE (e.g., `N`, `PRECISION`, `WINDOW_SIZE`)
- **Functions:** snake_case (e.g., `read_input()`, `getBasePairList()`)
- **Local variables:** lowercase (e.g., `i`, `j`, `energy`, `sequence`)
- **Constants:** UPPERCASE with underscore (e.g., `MAX_INTERIOR_DIST`)
- **Types/structs:** Suffix with `T` (e.g., `paramT`, `model_detailsT`)

### Memory Management
- Use `xcalloc()` wrapper for safe allocation
- Manual memory management with `free()` (no RAII/smart pointers)
- Dynamic 2D arrays frequently used

### Data Representations
- **Base pair list:** Array where `bpList[i]` = position paired with i (-1 if unpaired)
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
2. **Symmetry Optimization:** Only computes ceil(runLength/2) roots due to complex conjugate symmetry
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
1. Maintain C++98 compatibility - avoid C++11+ features
2. Use the established naming conventions
3. Preserve manual memory management patterns
4. Keep energy calculation precision in mind
5. Test with both Turner 1999 and 2004 parameters

### When Adding Features
1. Add new functions to appropriate files (delta.cpp for algorithm, misc.cpp for utilities)
2. Declare new functions in corresponding header files
3. Update CMakeLists.txt if adding new source files
4. Ensure OpenMP compatibility for parallelization

### Common Tasks
- **Build:** `cd build && cmake .. && make`
- **Run tests:** No formal test suite; verify with sample RNA sequences
- **Clean build:** `cd build && rm -rf *`

### Energy Parameters
- Default: `rna_turner2004.par`
- Alternative: `rna_turner1999.par`
- Format: ViennaRNA 2.x parameter format
- Installed to `bin/` directory alongside executable

### Important Notes
- The code uses Variable Length Arrays (VLAs) which are a C99/GCC extension
- OpenMP is required for parallel execution
- FFTW library handles the core FFT computations
- ViennaRNA library provides energy parameter handling functions
