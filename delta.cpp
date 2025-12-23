#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fftw3.h>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <memory>
#include "delta.h"
#include "misc.h"
#include "energy_par.h"
#include "params.h"
#include "parameter_parser.h"

// Constants
constexpr int STRUCTURE_COUNT = 1;
constexpr int MIN_PAIR_DIST = 3;
constexpr int MAX_INTERIOR_DIST = 30;
constexpr int FFTW_REAL = 0;
constexpr int FFTW_IMAG = 1;
constexpr bool FFTBOR_DEBUG = false;

// Complex number constants
const dcomplex ZERO_C(0.0, 0.0);
const dcomplex ONE_C(1.0, 0.0);

// Function-like macros (kept as macros due to external variable dependencies)
#define WINDOW_SIZE(i) (MIN_WINDOW_SIZE + i)
#define NUM_WINDOWS (WINDOW_SIZE - MIN_WINDOW_SIZE + 1)
#define ROOT_POW(i, pow, n) (rootsOfUnity[(i * pow) % (n + 1)])
#define ENERGY_DEBUG (0 && !root)

extern int    PF, N, PRECISION, WINDOW_SIZE, MIN_WINDOW_SIZE;
extern double temperature;
extern char   *ENERGY;
extern fftbor::ParamPtr P;

void neighbours(const char* inputSequence, const int* bpList) {
  int i, root, runLength = 0;
  const int sequenceLength = strlen(inputSequence);
  const double RT = 0.0019872370936902486 * (temperature + 273.15) * 100;

  const char* energyfile = ENERGY;

  // Use smart pointers and vectors for automatic memory management
  auto sequence = std::make_unique<char[]>(sequenceLength + 2);
  std::vector<int> intSequence(sequenceLength + 1);

  // Load energy parameters
  read_parameter_file(energyfile);
  P = scale_parameters();
  P->model_details.special_hp = 1;

  // Make necessary versions of the sequence for efficient processing
  sequence[0] = '@';
  strncpy(sequence.get() + 1, inputSequence, sequenceLength);
  sequence[sequenceLength + 1] = '\0';
  translateToIntSequence(inputSequence, intSequence.data());

  // Initialize canBasePair lookup for all 4x4 nt. combinations
  auto canBasePair = fftbor::make_int_matrix_2d(5, 5);
  initializeCanBasePair(canBasePair);

  // Initialize numBasePairs lookup for all (i, j) combinations
  auto numBasePairs = fftbor::make_int_matrix_2d(sequenceLength + 1, sequenceLength + 1);
  initializeBasePairCounts(numBasePairs, bpList, sequenceLength);

  // Determine max bp distance to determine how many roots of unity to generate
  for (i = 1; i <= sequenceLength; ++i) {
    runLength += (bpList[i] > i ? 1 : 0);
  }
  runLength += floor((sequenceLength - MIN_PAIR_DIST) / 2);

  // Use vectors for 2D matrices
  auto Z = fftbor::make_complex_matrix_2d(sequenceLength + 1, sequenceLength + 1);
  auto ZB = fftbor::make_complex_matrix_2d(sequenceLength + 1, sequenceLength + 1);
  auto ZM = fftbor::make_complex_matrix_2d(sequenceLength + 1, sequenceLength + 1);
  auto ZM1 = fftbor::make_complex_matrix_2d(sequenceLength + 1, sequenceLength + 1);

  // 3D solutions array
  fftbor::ComplexMatrix3D solutions(NUM_WINDOWS);
  for (i = 0; i < NUM_WINDOWS; ++i) {
    solutions[i].resize(sequenceLength - WINDOW_SIZE(i) + 2);
    for (int j = 1; j <= sequenceLength - WINDOW_SIZE(i) + 1; ++j) {
      solutions[i][j].resize(runLength + 1, ZERO_C);
    }
  }

  std::vector<dcomplex> rootsOfUnity(runLength + 1);
  populateMatrices(Z, ZB, ZM, ZM1, solutions, rootsOfUnity, sequenceLength, runLength);

  // Start main recursions (root <= ceil(runLength / 2.0) is an optimization for roots of unity)
  for (root = 0; root <= ceil(runLength / 2.0); ++root) {
    evaluateZ(root, Z, ZB, ZM, ZM1, solutions, rootsOfUnity,
              inputSequence, sequence.get(), intSequence.data(), bpList,
              canBasePair, numBasePairs, sequenceLength, runLength, RT);
  }

  if (FFTBOR_DEBUG) {
    std::cout << std::endl;
    printf("Number of structures: %.0f\n", Z[sequenceLength][1].real());
  }

  // Convert point-value solutions to coefficient form w/ inverse DFT
  populateRemainingRoots(solutions, sequenceLength, runLength, root);
  solveSystem(solutions, sequence.get(), bpList, sequenceLength, runLength);

  // Smart pointers and vectors clean up automatically - no manual cleanup needed!
  // P will be reset when neighbours is called again or on program exit
}

void evaluateZ(int root,
               fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
               fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
               fftbor::ComplexMatrix3D& solutions, const std::vector<dcomplex>& rootsOfUnity,
               const char* inputSequence, const char* sequence,
               const int* intSequence, const int* bpList,
               const fftbor::IntMatrix2D& canBasePair, const fftbor::IntMatrix2D& numBasePairs,
               int sequenceLength, int runLength, double RT) {
  int i, j, k, l, d, delta;
  double energy;

  flushMatrices(Z, ZB, ZM, ZM1, sequenceLength);

  if (ENERGY_DEBUG) {
    printf("RT: %f\n", RT / 100);
  }

  for (d = MIN_PAIR_DIST + 1; d < sequenceLength; ++d) {
    for (i = 1; i <= sequenceLength - d; ++i) {
      j = i + d;

      if (canBasePair[intSequence[i]][intSequence[j]]) {
        // ****************************************************************************
        // Solve ZB
        // ****************************************************************************
        // In a hairpin, [i + 1, j - 1] unpaired.
        // Pass pointer to position i-1 (0-indexed) so hairpinloop can read the correct loop motif
        energy    = hairpinloop(i, j, canBasePair[intSequence[i]][intSequence[j]], intSequence[i + 1], intSequence[j - 1], inputSequence + i - 1);
        delta     = numBasePairs[i][j] + jPairedTo(i, j, bpList);
        ZB[i][j] += ROOT_POW(root, delta, runLength) * exp(-energy / RT);

        if (ENERGY_DEBUG) {
          printf("%+f: GetHairpinEnergy(%c (%d), %c (%d));\n", energy / 100, sequence[i], i, sequence[j], j);
        }

        if (STRUCTURE_COUNT) {
          ZB[j][i] += 1;
        }

        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < min2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          for (l = max2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) {
            if (canBasePair[intSequence[k]][intSequence[l]]) {
              energy    = interiorloop(i, j, k, l, canBasePair[intSequence[i]][intSequence[j]], canBasePair[intSequence[l]][intSequence[k]], intSequence[i + 1], intSequence[l + 1], intSequence[j - 1], intSequence[k - 1]);
              delta     = numBasePairs[i][j] - numBasePairs[k][l] + jPairedTo(i, j, bpList);
              ZB[i][j] += ZB[k][l] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

              if (ENERGY_DEBUG) {
                printf("%+f: GetInteriorStackingAndBulgeEnergy(%c (%d), %c (%d), %c (%d), %c (%d));\n", energy / 100, sequence[i], i, sequence[j], j, sequence[k], k, sequence[l], l);
              }

              if (STRUCTURE_COUNT) {
                ZB[j][i] += ZB[l][k];
              }
            }
          }
        }

        for (k = i + MIN_PAIR_DIST + 3; k < j - MIN_PAIR_DIST - 1; ++k) {
          energy    = P->MLclosing + P->MLintern[canBasePair[intSequence[i]][intSequence[j]]];
          delta     = numBasePairs[i][j] - numBasePairs[i + 1][k - 1] - numBasePairs[k][j - 1] + jPairedTo(i, j, bpList);
          ZB[i][j] += ZM[i + 1][k - 1] * ZM1[k][j - 1] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

          if (ENERGY_DEBUG) {
            printf("%+f: MultiloopA + MultiloopB;\n", energy / 100);
          }

          if (STRUCTURE_COUNT) {
            ZB[j][i] += ZM[k - 1][i + 1] * ZM1[j - 1][k];
          }
        }
      }

      // ****************************************************************************
      // Solve ZM1
      // ****************************************************************************
      for (k = i + MIN_PAIR_DIST + 1; k <= j; ++k) {
        if (canBasePair[intSequence[i]][intSequence[k]]) {
          energy     = P->MLbase * (j - k) + P->MLintern[canBasePair[intSequence[i]][intSequence[k]]];
          delta      = numBasePairs[i][j] - numBasePairs[i][k];
          ZM1[i][j] += ZB[i][k] * exp(-energy / RT);

          if (ENERGY_DEBUG) {
            printf("%+f: MultiloopB + MultiloopC * (%d - %d);\n", energy / 100, j, k);
          }

          if (STRUCTURE_COUNT) {
            ZM1[j][i] += ZB[k][i];
          }
        }
      }

      // ****************************************************************************
      // Solve ZM
      // ****************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        energy    = P->MLbase * (k - i);
        delta     = numBasePairs[i][j] - numBasePairs[k][j];
        ZM[i][j] += ZM1[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

        if (ENERGY_DEBUG) {
          printf("%+f: MultiloopC * (%d - %d);\n", energy / 100, k, i);
        }

        if (STRUCTURE_COUNT) {
          ZM[j][i] += ZM1[j][k];
        }

        if (k > i + MIN_PAIR_DIST + 1) {
          energy    = P->MLintern[canBasePair[intSequence[k]][intSequence[j]]];
          delta     = numBasePairs[i][j] - numBasePairs[i][k - 1] - numBasePairs[k][j];
          ZM[i][j] += ZM[i][k - 1] * ZM1[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

          if (ENERGY_DEBUG) {
            printf("%+f: MultiloopB;\n", energy / 100);
          }

          if (STRUCTURE_COUNT) {
            ZM[j][i] += ZM[k - 1][i] * ZM1[j][k];
          }
        }
      }

      // **************************************************************************
      // Solve Z
      // **************************************************************************
      delta    = jPairedIn(i, j, bpList);
      Z[i][j] += Z[i][j - 1] * ROOT_POW(root, delta, runLength);

      if (STRUCTURE_COUNT) {
        Z[j][i] += Z[j - 1][i];
      }

      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        if (canBasePair[intSequence[k]][intSequence[j]]) {
          energy = canBasePair[intSequence[k]][intSequence[j]] > 2 ? P->TerminalAU : 0;

          if (ENERGY_DEBUG) {
            printf("%+f: %c-%c == (2 || 3) ? 0 : GUAU_penalty;\n", energy / 100, sequence[k], sequence[j]);
          }

          if (k == i) {
            delta    = numBasePairs[i][j] - numBasePairs[k][j];
            Z[i][j] += ZB[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

            if (STRUCTURE_COUNT) {
              Z[j][i] += ZB[j][k];
            }
          } else {
            delta    = numBasePairs[i][j] - numBasePairs[i][k - 1] - numBasePairs[k][j];
            Z[i][j] += Z[i][k - 1] * ZB[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

            if (STRUCTURE_COUNT) {
              Z[j][i] += Z[k - 1][i] * ZB[j][k];
            }
          }
        }
      }
    }
  }

  for (i = 0; i < NUM_WINDOWS; ++i) {
    for (j = 1; j <= sequenceLength - WINDOW_SIZE(i) + 1; ++j) {
      solutions[i][j][root] = Z[j][j + WINDOW_SIZE(i) - 1];
    }
  }

  if (FFTBOR_DEBUG) {
    std::cout << "." << std::flush;
  }
}

void solveSystem(fftbor::ComplexMatrix3D& solutions, const char* sequence, const int* structure, int sequenceLength, int runLength) {
  char precisionFormat[20];
  int i, j, k;
  double scalingFactor, sum;

  snprintf(precisionFormat, sizeof(precisionFormat), "%%d\t%%.0%df\n", PRECISION ? PRECISION : std::numeric_limits<double>::digits10);

  // Use FFTW's allocation for proper alignment (fftw_complex is double[2], can't use std::vector)
  struct FftwDeleter {
    void operator()(fftw_complex* p) const { fftw_free(p); }
  };
  std::unique_ptr<fftw_complex[], FftwDeleter> signal(fftw_alloc_complex(runLength + 1));
  std::unique_ptr<fftw_complex[], FftwDeleter> result(fftw_alloc_complex(runLength + 1));

  if (!signal || !result) {
    fprintf(stderr, "Error: Failed to allocate FFTW arrays (size: %d)\n", runLength + 1);
    exit(1);
  }

  fftw_plan plan = fftw_plan_dft_1d(runLength + 1, signal.get(), result.get(), FFTW_FORWARD, FFTW_ESTIMATE);

  if (!plan) {
    fprintf(stderr, "Error: FFTW plan creation failed for size %d\n", runLength + 1);
    exit(1);
  }

  for (i = 0; i < NUM_WINDOWS; ++i) {
    for (j = 1; j <= sequenceLength - WINDOW_SIZE(i) + 1; ++j) {
      sum           = 0;
      scalingFactor = solutions[i][j][0].real();

      if (!(MIN_WINDOW_SIZE == WINDOW_SIZE && WINDOW_SIZE == N)) {
        printf("Window size:           %d\n", WINDOW_SIZE(i));
        printf("Window starting index: %d\n", j);

        printf("Sequence  (%d, %d): ", j, j + WINDOW_SIZE(i) - 1);
        for (k = j; k <= j + WINDOW_SIZE(i) - 1; k++) {
          printf("%c", sequence[k]);
        }
        printf("\n");

        printf("Structure (%d, %d): ", j, j + WINDOW_SIZE(i) - 1);
        for (k = j; k <= j + WINDOW_SIZE(i) - 1; k++) {
          printf("%c", structure[k] < 0 ? '.' : (structure[k] > k ? '(' : ')'));
        }
        printf("\n");
      }

      for (k = 0; k <= runLength; k++) {
        signal[k][FFTW_REAL] = (pow(10, PRECISION) * solutions[i][j][k].real()) / scalingFactor;
        signal[k][FFTW_IMAG] = (pow(10, PRECISION) * solutions[i][j][k].imag()) / scalingFactor;
      }

      fftw_execute(plan);

      printf("k\tp(k)\n");

      for (k = 0; k <= runLength; k++) {
        if (PRECISION == 0) {
          solutions[i][j][k] = dcomplex(result[k][FFTW_REAL] / (runLength + 1), 0);
        } else {
          solutions[i][j][k] = dcomplex(pow(10.0, -PRECISION) * static_cast<int>(result[k][FFTW_REAL] / (runLength + 1)), 0);
        }

        sum += solutions[i][j][k].real();

        printf(precisionFormat, k, solutions[i][j][k].real());
      }

      if (FFTBOR_DEBUG) {
        printf("Scaling factor (Z{%d, %d}): %.15f\n", j, j + WINDOW_SIZE(i) - 1, scalingFactor);
        std::cout << "Sum: " << sum << std::endl << std::endl;
      }
    }
  }

  fftw_destroy_plan(plan);
}

int jPairedTo(int i, int j, const int* basePairs) {
  return basePairs[i] == j ? -1 : 1;
}

int jPairedIn(int i, int j, const int* basePairs) {
  return basePairs[j] >= i && basePairs[j] < j ? 1 : 0;
}

void populateRemainingRoots(fftbor::ComplexMatrix3D& solutions, int sequenceLength, int runLength, int lastRoot) {
  int i, j, k, root;

  for (i = 0; i < NUM_WINDOWS; ++i) {
    for (j = 1; j <= sequenceLength - WINDOW_SIZE(i) + 1; ++j) {
      root = lastRoot;

      if (runLength % 2) {
        k = root - 2;
      } else {
        k = root - 1;
      }

      for (; root <= runLength && k > 0; --k, ++root) {
        solutions[i][j][root] = dcomplex(solutions[i][j][k].real(), -solutions[i][j][k].imag());
      }
    }
  }
}

void populateMatrices(fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                      fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                      fftbor::ComplexMatrix3D& solutions, std::vector<dcomplex>& rootsOfUnity,
                      int sequenceLength, int runLength) {
  // Matrices are already sized correctly - just initialize rootsOfUnity
  for (int i = 0; i <= runLength; ++i) {
    rootsOfUnity[i] = dcomplex(cos(2 * M_PI * i / (runLength + 1)), sin(2 * M_PI * i / (runLength + 1)));
  }
}

void flushMatrices(fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                   fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                   int sequenceLength) {
  for (int i = 0; i <= sequenceLength; ++i) {
    for (int j = 0; j <= sequenceLength; ++j) {
      if (i > 0 && j > 0 && abs(j - i) <= MIN_PAIR_DIST) {
        Z[i][j] = ONE_C;
      } else {
        Z[i][j] = ZERO_C;
      }

      ZB[i][j]  = ZERO_C;
      ZM[i][j]  = ZERO_C;
      ZM1[i][j] = ZERO_C;
    }
  }
}

int numbp(int i, int j, const int* bpList) {
  int n = 0;
  for (int k = i; k <= j; k++) {
    if (k < bpList[k] && bpList[k] <= j) {
      n++;
    }
  }
  return n;
}

int bn(char A) {
  if (A == 'A') return 1;
  if (A == 'C') return 2;
  if (A == 'G') return 3;
  if (A == 'U') return 4;
  return 0;
}

void initializeCanBasePair(fftbor::IntMatrix2D& canBasePair) {
  // A = 1, C = 2, G = 3, U = 4
  canBasePair[1][4] = 5;
  canBasePair[4][1] = 6;
  canBasePair[2][3] = 1;
  canBasePair[3][2] = 2;
  canBasePair[3][4] = 3;
  canBasePair[4][3] = 4;
}

void translateToIntSequence(const char* a, int* intSequence) {
  int len = strlen(a);
  intSequence[0] = len;
  for (int i = 0; i < len; i++) {
    intSequence[i + 1] = bn(a[i]);
  }
}

void initializeBasePairCounts(fftbor::IntMatrix2D& numBasePairs, const int* bpList, int n) {
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      numBasePairs[i][j] = 0;
    }
  }
  for (int d = MIN_PAIR_DIST + 1; d < n; d++) {
    for (int i = 1; i <= n - d; i++) {
      int j = i + d;
      numBasePairs[i][j] = numbp(i, j, bpList);
    }
  }
}

// Calculate hairpin loop energy.
// Parameters:
//   i, j: 1-indexed positions of the closing base pair
//   type: base pair type (1-6 for valid pairs)
//   si1, sj1: integer encodings of bases at i+1 and j-1
//   string: pointer to the START of the hairpin loop in the sequence (position i-1 in 0-indexed)
//           The caller must ensure string points to at least (size + 2) valid characters.
double hairpinloop(int i, int j, int type, short si1, short sj1, const char* string) {
  double energy;
  int size = j - i - 1;

  energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30] + (int)(P->lxc * log((size) / 30.));
  if (P->model_details.special_hp) {
    if (size == 4) {
      // Tetraloop: 6 characters (closing pair + 4 loop bases)
      char tl[7] = {0};
      strncpy(tl, string, 6);
      tl[6] = '\0';  // Ensure null termination
      char* ts = strstr(P->Tetraloops, tl);
      if (ts) {
        return P->Tetraloop_E[(ts - P->Tetraloops) / 7];
      }
    } else if (size == 6) {
      // Hexaloop: 8 characters (closing pair + 6 loop bases)
      char tl[9] = {0};
      strncpy(tl, string, 8);
      tl[8] = '\0';  // Ensure null termination
      char* ts = strstr(P->Hexaloops, tl);
      if (ts) {
        return P->Hexaloop_E[(ts - P->Hexaloops) / 9];
      }
    } else if (size == 3) {
      // Triloop: 5 characters (closing pair + 3 loop bases)
      char tl[6] = {0};
      strncpy(tl, string, 5);
      tl[5] = '\0';  // Ensure null termination
      char* ts = strstr(P->Triloops, tl);
      if (ts) {
        return P->Triloop_E[(ts - P->Triloops) / 6];
      }
      return energy + (type > 2 ? P->TerminalAU : 0);
    }
  }
  energy += P->mismatchH[type][si1][sj1];

  return energy;
}

double interiorloop(int i, int j, int k, int l, int type, int type_2, short si1, short sq1, short sj1, short sp1) {
  int nl, ns;
  double energy;
  int n1 = k - i - 1;
  int n2 = j - l - 1;
  energy = INF;

  if (n1 > n2) {
    nl = n1;
    ns = n2;
  } else {
    nl = n2;
    ns = n1;
  }

  if (nl == 0) {
    return P->stack[type][type_2];
  }

  if (ns == 0) {
    energy = (nl <= MAXLOOP) ? P->bulge[nl] : (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
    if (nl == 1) {
      energy += P->stack[type][type_2];
    } else {
      if (type > 2) energy += P->TerminalAU;
      if (type_2 > 2) energy += P->TerminalAU;
    }
    return energy;
  } else {
    if (ns == 1) {
      if (nl == 1) {
        return P->int11[type][type_2][si1][sj1];
      }
      if (nl == 2) {
        if (n1 == 1) {
          energy = P->int21[type][type_2][si1][sq1][sj1];
        } else {
          energy = P->int21[type_2][type][sq1][si1][sp1];
        }
        return energy;
      } else {
        energy = (nl + 1 <= MAXLOOP) ? P->internal_loop[nl + 1] : (P->internal_loop[30] + (int)(P->lxc * log((nl + 1) / 30.)));
        energy += std::min(MAX_NINIO, (nl - ns) * P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    } else if (ns == 2) {
      if (nl == 2) {
        return P->int22[type][type_2][si1][sp1][sq1][sj1];
      } else if (nl == 3) {
        energy = P->internal_loop[5] + P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }
    }
    {
      energy = (n1 + n2 <= MAXLOOP) ? P->internal_loop[n1 + n2] : (P->internal_loop[30] + (int)(P->lxc * log((n1 + n2) / 30.)));
      energy += std::min(MAX_NINIO, (nl - ns) * P->ninio[2]);
      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}
