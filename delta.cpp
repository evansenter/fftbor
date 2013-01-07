#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include "delta.h"
#include "misc.h"
#include <fftw3.h>
#include "energy_const.h"
#include "energy_par.h"
#include <iostream>
#include <limits>

#define STRUCTURE_COUNT 1
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define FFTW_REAL 0
#define FFTW_IMAG 1
#define DELTA2D(expression1, expression2, n) (((expression1) * (n)) + (expression2))
#define ROOT_POW(i, pow, n) (rootsOfUnity[(i * pow) % n])
#define FFTBOR_DEBUG 1
#define ENERGY_DEBUG (0 && !root)
#define TABLE_HEADERS 1

extern int    PF, N, PRECISION, WINDOW_SIZE, MIN_WINDOW_SIZE;
extern double temperature;
extern char   *ENERGY;
extern paramT *P;

extern "C" void read_parameter_file(const char energyfile[]);

void neighbors(char *inputSequence, int **bpList) {
  int i, root, rowLength, runLength, sequenceLength = strlen(inputSequence), inputStructureDist = 0;
  double RT = 0.0019872370936902486 * (temperature + 273.15) * 100; // 0.01 * (kcal K) / mol

  char *energyfile    = ENERGY;
  char *sequence      = new char[sequenceLength + 1];
  int  *intSequence   = (int *)xcalloc(sequenceLength + 1, sizeof(int));
  int ***numBasePairs = new int**[2];
  int **canBasePair;
  
  // Load energy parameters.
  read_parameter_file(energyfile);
  P = scale_parameters();

  // Make necessary versions of the sequence for efficient processing.
  sequence[0] = '@';
  strncpy(sequence + 1, inputSequence, sequenceLength);
  translateToIntSequence(inputSequence, intSequence);
  
  // Initialize canBasePair lookup for all 4x4 nt. combinations.
  canBasePair = (int **)xcalloc(5, sizeof(int *));
  for (i = 0; i < 5; ++i) {
    canBasePair[i] = (int *)xcalloc(5, sizeof(int));
  }
  initializeCanBasePair(canBasePair);
  
  // Initialize numBasePairs lookup for all (i, j) combinations.
  numBasePairs[0] = (int **)xcalloc(sequenceLength + 1, sizeof(int *));
  numBasePairs[1] = (int **)xcalloc(sequenceLength + 1, sizeof(int *));
  for (i = 1; i <= sequenceLength; i++) {
    numBasePairs[0][i] = (int *)xcalloc(sequenceLength + 1, sizeof(int));
    numBasePairs[1][i] = (int *)xcalloc(sequenceLength + 1, sizeof(int));
  }
  initializeBasePairCounts(numBasePairs[0], bpList[0], sequenceLength);
  initializeBasePairCounts(numBasePairs[1], bpList[1], sequenceLength);
  
  // Initialize the inputStructureDist value as bp_dist(S_a, S_b)
  for (i = 1; i <= sequenceLength; ++i) {
    inputStructureDist += (bpList[0][i] > i && bpList[0][i] != bpList[1][i] ? 1 : 0);
    inputStructureDist += (bpList[1][i] > i && bpList[1][i] != bpList[0][i] ? 1 : 0);
  }
  
  // Note: rowLength = (least even number >= sequenceLength) + 1 (for a seq. of length 9 rowLength = (0..10).length = 11)
  rowLength = (sequenceLength % 2 ? sequenceLength + 1 : sequenceLength) + 1;
  // Note: runLength = least number div. 4 >= (rowLength ^ 2 + 1) (for a seq. of length 9 runLength = (11 ^ 2 + 1) + 2 = 124)
  runLength = ((pow(rowLength, 2) + 1) + (int)(pow(rowLength, 2) + 1) % 4);
  
  dcomplex **Z            = new dcomplex*[sequenceLength + 1];
  dcomplex **ZB           = new dcomplex*[sequenceLength + 1];
  dcomplex **ZM           = new dcomplex*[sequenceLength + 1];
  dcomplex **ZM1          = new dcomplex*[sequenceLength + 1];
  dcomplex *solutions     = new dcomplex[runLength];
  dcomplex *rootsOfUnity  = new dcomplex[runLength];
  
  populateMatrices(Z, ZB, ZM, ZM1, solutions, rootsOfUnity, sequenceLength, runLength);
  
  if (FFTBOR_DEBUG) {
    printf("rowLength:          %d\n", rowLength);
    printf("runLength:          %d\n", runLength);
    printf("sequenceLength:     %d\n", sequenceLength);
    printf("inputStructureDist: %d\n", inputStructureDist);
    printf("Roots of unity:\n");
    for (root = 0; root < runLength; ++root) {
      printf("%d: %+f %+fi\n", root, rootsOfUnity[root].real(), rootsOfUnity[root].imag());
    }
  }
  
  // Start main recursions.
  for (root = 0; root <= runLength / 4; ++root) {
    evaluateZ(2 * root, Z, ZB, ZM, ZM1, solutions, rootsOfUnity, inputSequence, sequence, intSequence, bpList, canBasePair, numBasePairs, inputStructureDist, sequenceLength, rowLength, runLength, RT);
  }
  
  if (FFTBOR_DEBUG) {
    std::cout << std::endl;
    printf("Number of structures: %.0f\n", Z[sequenceLength][1].real());
  }

  // Convert point-value solutions to coefficient form w/ inverse DFT.
  populateRemainingRoots(solutions, sequenceLength, runLength);
  
  if (inputStructureDist % 2) {
    // Odd case.
    for (root = 0; root < runLength; ++root) {      
      solutions[root] = rootsOfUnity[(runLength - root) % runLength] * solutions[root];
    }
  }
  
  solveSystem(solutions, sequence, bpList, sequenceLength, rowLength, runLength, inputStructureDist);
  
  free(intSequence);
}

void evaluateZ(int root, dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, dcomplex *solutions, dcomplex *rootsOfUnity, char *inputSequence, char *sequence, int *intSequence, int *bpList[2], int *canBasePair[5], int **numBasePairs[2], int inputStructureDist, int sequenceLength, int rowLength, int runLength, double RT) {
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
        energy = hairpinloop(i, j, canBasePair[intSequence[i]][intSequence[j]], intSequence, inputSequence);
        delta  = DELTA2D(
          numBasePairs[0][i][j] + jPairedTo(i, j, bpList[0]), 
          numBasePairs[1][i][j] + jPairedTo(i, j, bpList[1]), 
          rowLength
        );
        ZB[i][j] += ROOT_POW(root, delta, runLength) * exp(-energy / RT);

        if (ENERGY_DEBUG) {
          printf("%+f: GetHairpinEnergy(%c (%d), %c (%d)); Delta = %d\n", energy / 100, sequence[i], i, sequence[j], j, delta);
        }

        if (STRUCTURE_COUNT) {
          ZB[j][i] += 1;
        }

        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < min2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          // There can't be more than 30 unpaired bases between i and k, and there must be room between k and j for l
          for (l = max2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) { 
            // l needs to at least have room to pair with k, and there can be at most 30 unpaired bases between (i, k) + (l, j), with l < j
            if (canBasePair[intSequence[k]][intSequence[l]]) {
              // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
              energy = interiorloop(i, j, k, l, canBasePair[intSequence[i]][intSequence[j]], canBasePair[intSequence[l]][intSequence[k]], intSequence);
              delta  = DELTA2D(
                numBasePairs[0][i][j] - numBasePairs[0][k][l] + jPairedTo(i, j, bpList[0]), 
                numBasePairs[1][i][j] - numBasePairs[1][k][l] + jPairedTo(i, j, bpList[1]), 
                rowLength
              );
              ZB[i][j] += ZB[k][l] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);
                
              if (ENERGY_DEBUG) {
                printf("%+f: GetInteriorStackingAndBulgeEnergy(%c (%d), %c (%d), %c (%d), %c (%d)); Delta = %d\n", energy / 100, sequence[i], i, sequence[j], j, sequence[k], k, sequence[l], l, delta);
              }

              if (STRUCTURE_COUNT) {
                ZB[j][i] += ZB[l][k];
              }
            }
          }
        }
        
        for (k = i + MIN_PAIR_DIST + 3; k < j - MIN_PAIR_DIST - 1; ++k) {
          // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1)
          energy = P->MLclosing + P->MLintern[canBasePair[intSequence[i]][intSequence[j]]];
          delta  = DELTA2D(
            numBasePairs[0][i][j] - numBasePairs[0][i + 1][k - 1] - numBasePairs[0][k][j - 1] + jPairedTo(i, j, bpList[0]), 
            numBasePairs[1][i][j] - numBasePairs[1][i + 1][k - 1] - numBasePairs[1][k][j - 1] + jPairedTo(i, j, bpList[1]), 
            rowLength
          );
          ZB[i][j] += ZM[i + 1][k - 1] * ZM1[k][j - 1] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);
                 
          if (ENERGY_DEBUG) {
            printf("%+f: MultiloopA + MultiloopB; Delta = %d\n", energy / 100, delta);
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
        // k is the closing base pairing with i of a single component within the range [i, j]
        if (canBasePair[intSequence[i]][intSequence[k]]) {
          energy = P->MLbase * (j - k) + P->MLintern[canBasePair[intSequence[i]][intSequence[k]]];
          delta  = DELTA2D(
            numBasePairs[0][i][j] - numBasePairs[0][i][k],
            numBasePairs[1][i][j] - numBasePairs[1][i][k],
            rowLength
          );
          
          ZM1[i][j] += ZB[i][k] * exp(-energy / RT);
          
          if (ENERGY_DEBUG) {
            printf("%+f: MultiloopB + MultiloopC * (%d - %d); Delta = %d\n", energy / 100, j, k, delta);
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
        // Only one stem.
        energy = P->MLbase * (k - i);
        delta  = DELTA2D(
          numBasePairs[0][i][j] - numBasePairs[0][k][j],
          numBasePairs[1][i][j] - numBasePairs[1][k][j],
          rowLength
        );
        ZM[i][j] += ZM1[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

        if (ENERGY_DEBUG) {
          printf("%+f: MultiloopC * (%d - %d); Delta = %d\n", energy / 100, k, i, delta);
        }

        if (STRUCTURE_COUNT) {
          ZM[j][i] += ZM1[j][k];
        }

        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
          energy = P->MLintern[canBasePair[intSequence[k]][intSequence[j]]];
          delta  = DELTA2D(
            numBasePairs[0][i][j] - numBasePairs[0][i][k - 1] - numBasePairs[0][k][j],
            numBasePairs[1][i][j] - numBasePairs[1][i][k - 1] - numBasePairs[1][k][j],
            rowLength
          );
          ZM[i][j] += ZM[i][k - 1] * ZM1[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

          if (ENERGY_DEBUG) {
            printf("%+f: MultiloopB; Delta = %d\n", energy / 100, delta);
          }

          if (STRUCTURE_COUNT) {
            ZM[j][i] += ZM[k - 1][i] * ZM1[j][k];
          }
        }
      }
        
      // **************************************************************************
      // Solve Z
      // **************************************************************************
      delta = DELTA2D(
        jPairedIn(i, j, bpList[0]),
        jPairedIn(i, j, bpList[1]),
        rowLength
      );
      Z[i][j] += Z[i][j - 1] * ROOT_POW(root, delta, runLength);

      if (STRUCTURE_COUNT) {
        Z[j][i] += Z[j - 1][i];
      }

      for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
        // (k, j) is the rightmost base pair in (i, j)
        if (canBasePair[intSequence[k]][intSequence[j]]) {
          energy = canBasePair[intSequence[k]][intSequence[j]] > 2 ? TerminalAU : 0;

          if (ENERGY_DEBUG) {
            printf("%+f: %c-%c == (2 || 3) ? 0 : GUAU_penalty; Delta = %d\n", energy / 100, sequence[k], sequence[j], delta);
          }
            
          if (k == i) {
            delta = DELTA2D(
              numBasePairs[0][i][j] - numBasePairs[0][k][j],
              numBasePairs[1][i][j] - numBasePairs[1][k][j],
              rowLength
            );
            Z[i][j] += ZB[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

            if (STRUCTURE_COUNT) {
              Z[j][i] += ZB[j][k];
            }
          } else {
            delta = DELTA2D(
              numBasePairs[0][i][j] - numBasePairs[0][i][k - 1] - numBasePairs[0][k][j],
              numBasePairs[1][i][j] - numBasePairs[1][i][k - 1] - numBasePairs[1][k][j],
              rowLength
            );
            Z[i][j] += Z[i][k - 1] * ZB[k][j] * ROOT_POW(root, delta, runLength) * exp(-energy / RT);

            if (STRUCTURE_COUNT) {
              Z[j][i] += Z[k - 1][i] * ZB[j][k];
            }
          }
        }
      }
    }
  }
  
  solutions[root] = Z[1][sequenceLength];

  if (FFTBOR_DEBUG) {
    std::cout << root << ' ' << std::flush;
  }
}

void solveSystem(dcomplex *solutions, char *sequence, int **structure, int sequenceLength, int rowLength, int runLength, int inputStructureDist) {
  char precisionFormat[20];
  char header[5 * rowLength]; // This is enough space for sequences up to length 999
  int i, solutionLength = pow(rowLength, 2);
  double scalingFactor, sum;
  double *probabilities = (double *)xcalloc(solutionLength, sizeof(double));
  
  if (TABLE_HEADERS) {
    sprintf(precisionFormat, "\t%%+.0%df", PRECISION ? PRECISION : std::numeric_limits<double>::digits10);
  } else {
    sprintf(precisionFormat, "%%+.0%df\t", PRECISION ? PRECISION : std::numeric_limits<double>::digits10);
  }
  
  for (i = 0; i < rowLength; ++i) {
    sprintf(&header[4 * i], "\t%-3d", i);
  }
  
  if (FFTBOR_DEBUG) {
    printf("Solutions (after populating remaining roots):\n");
    for (i = 0; i < runLength; ++i) {
      printf("%d: %+f %+fi\n", i, solutions[i].real(), solutions[i].imag());
    }
  }
  
  fftw_complex signal[runLength / 2];
  fftw_complex result[runLength / 2];
  
  fftw_plan plan = fftw_plan_dft_1d(runLength / 2, signal, result, FFTW_BACKWARD, FFTW_ESTIMATE);
  sum            = 0;
  scalingFactor  = solutions[0].real();
  
  // For some reason it's much more numerically stable to set the signal via real / imag components separately.
  for (i = 0; i < runLength / 2; ++i) {
    // Convert point-value solutions of Z(root) to 10^PRECISION * Z(root) / Z
    signal[i][FFTW_REAL] = (pow(10, PRECISION) * solutions[2 * i].real()) / scalingFactor;
    signal[i][FFTW_IMAG] = (pow(10, PRECISION) * solutions[2 * i].imag()) / scalingFactor;
  }
  
  // Calculate transform, coefficients are in fftw_complex result array.
  fftw_execute(plan);
  
  for (i = 0; i < runLength / 2; ++i) {
    // Truncate to user-specified precision, default is 4 and if set to 0, no truncation occurs (dangerous).
    if (PRECISION == 0) {
      solutions[i] = dcomplex(result[i][FFTW_REAL] / (runLength / 2), 0);
    } else {
      solutions[i] = dcomplex(pow(10.0, -PRECISION) * static_cast<int>(result[i][FFTW_REAL] / (runLength / 2)), 0);
    }
    
    if (inputStructureDist % 2) {
      // Odd case
      probabilities[2 * i + 1] = solutions[i].real();
    } else {
      // Even case
      probabilities[2 * i] = solutions[i].real();
    }
    
    sum += solutions[i].real();
  }
  
  if (TABLE_HEADERS) {
    std::cout << "Unspaced table:" << std::endl;
    std::cout << header << std::endl;
  }
  
  for (i = 0; i < solutionLength; ++i) {
    if (TABLE_HEADERS) {
      if (!i) {
        printf("0");
      } else if (!(i % rowLength)) {
        printf("\n%d", i / rowLength);
      }
    } else if (!(i % rowLength)) {
      printf("\n");
    }
    
    if (i < runLength / 2) {
      printf(precisionFormat, solutions[i].real());
    } else {
      printf(precisionFormat, 0);
    }
  }
      
  printf("\n\n");
  
  if (TABLE_HEADERS) {
    std::cout << "Spaced table:" << std::endl;
    std::cout << header << std::endl;
  }
  
  for (i = 0; i < solutionLength; ++i) {
    if (TABLE_HEADERS) {
      if (!i) {
        printf("0");
      } else if (!(i % rowLength)) {
        printf("\n%d", i / rowLength);
      }
    } else if (!(i % rowLength)) {
      printf("\n");
    }
    
    printf(precisionFormat, probabilities[i]);
  }
      
  printf("\n\n");
  
  if (FFTBOR_DEBUG) {
    for (i = 0; i < solutionLength; ++i) {
      // If the parity of (i + j) doesn't equal bp_dist(S_a, S_b) and p_{i, j} > 0, we have a problem (by the triangle inequality)
      if ((((i / rowLength) + (i % rowLength)) % 2) != (inputStructureDist % 2) && (int)(probabilities[i] * pow(10, PRECISION)) > 0) {
        printf("Warning: non-zero entry at (%d,\t%d):\t", i / rowLength, i % rowLength);
        printf(precisionFormat, probabilities[i]);
        printf("\n");
      }
    }
    
    printf("Scaling factor: %.15f\n", scalingFactor);
    std::cout << "Sum: " << sum << std::endl << std::endl;
  }
  
  fftw_destroy_plan(plan);
}

int jPairedTo(int i, int j, int *basePairs) {
  return basePairs[i] == j ? -1 : 1;
}

int jPairedIn(int i, int j, int *basePairs) {
  return basePairs[j] >= i && basePairs[j] < j ? 1 : 0;
}

void populateRemainingRoots(dcomplex *solutions, int sequenceLength, int runLength) {
  // Optimization leveraging complex conjugate of solutions to the polynomial.
  int i, l = runLength / 2 - 1, r = runLength / 2 + 1;
  
  if (FFTBOR_DEBUG) {
    printf("runLength: %d, runLength %% 2: %d\n\n", runLength, runLength % 2);
    
    printf("Solutions (before populating remaining roots):\n");
    for (i = 0; i < runLength; ++i) {
      printf("%d: %+f %+fi\n", i, solutions[i].real(), solutions[i].imag());
    }
  }
  
  for (i = 0; i < runLength / 2 && l - i >= 0 && r + i < runLength; ++i) {
    printf("i: %d, l: %d, r %d\n", i, l - i, r + i);
    solutions[r + i] = dcomplex(solutions[l - i].real(), -solutions[l - i].imag());
  }
}

void populateMatrices(dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, dcomplex *solutions, dcomplex *rootsOfUnity, int sequenceLength, int runLength) {
  int i;
  
  for (i = 0; i <= sequenceLength; ++i) {
    Z[i]   = new dcomplex[sequenceLength + 1];
    ZB[i]  = new dcomplex[sequenceLength + 1];
    ZM[i]  = new dcomplex[sequenceLength + 1];
    ZM1[i] = new dcomplex[sequenceLength + 1];
  }
  
  for (i = 0; i < runLength; ++i) {
    rootsOfUnity[i] = dcomplex(cos(-2 * M_PI * i / runLength), sin(-2 * M_PI * i / runLength));
  }
}

void flushMatrices(dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, int sequenceLength) {
  int i, j;
  
  for (i = 0; i <= sequenceLength; ++i) {
    for (j = 0; j <= sequenceLength; ++j) {
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

/* Number of base pairs in the region i to j in bpList */
int numbp(int i, int j, int *bpList) {
  int n=0;
  int k;
  for (k=i;k<=j;k++)
    if ( k<bpList[k] && bpList[k]<=j )
      n++;
  return n;
}

int bn(char A) {
/* Return the integer corresponding to the given base
   @ = 0  A = 1  C = 2  G = 3  U = 4 */
  if (A=='A')
    return 1;
  else if (A=='C')
    return 2;
  else if (A=='G')
    return 3;
  else if (A=='U')
    return 4;
  else
    return 0;
}

void initializeCanBasePair(int **canBasePair) {
  // A = 1, C = 2, G = 3, U = 4
  canBasePair[1][4] = 5;
  canBasePair[4][1] = 6;
  canBasePair[2][3] = 1;
  canBasePair[3][2] = 2;
  canBasePair[3][4] = 3;
  canBasePair[4][3] = 4;
}

void translateToIntSequence(char *a, int *intSequence) {
  int i;
  intSequence[0] = strlen(a);
  for (i=0;i<intSequence[0];i++)
    intSequence[i+1] = bn(a[i]);
}

void initializeBasePairCounts(int **numBasePairs, int *bpList, int n){
  int d,i,j;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      numBasePairs[i][j] = 0;
  for (d = MIN_PAIR_DIST+1; d < n; d++) {
    for (i = 1; i <= n - d; i++) {
      j = i + d;
      numBasePairs[i][j] = numbp(i,j,bpList);
    }
  }
}

double hairpinloop(int i, int j, int bp_type, int *intSequence, char *a) {
  // char *a is intended to be 0-indexed.
  double energy;
  energy = ((j-i-1) <= 30) ? P->hairpin[j-i-1] :
    P->hairpin[30]+(P->lxc*log((j-i-1)/30.));
  if ((j-i-1) >3) /* No mismatch for triloops */
    energy += P->mismatchH[bp_type]
      [intSequence[i+1]][intSequence[j-1]];
  else {
    char tl[6]={0}, *ts;
    strncpy(tl, a+i-1, 5);
    if ((ts=strstr(Triloops, tl))) 
      energy += P->Triloop_E[(ts - Triloops)/6];
    if (bp_type>2)
      energy += TerminalAU;
	}
  if ((j-i-1) == 4) { /* check for tetraloop bonus */
    char tl[7]={0}, *ts;
    strncpy(tl, a+i-1, 6);
    if ((ts=strstr(Tetraloops, tl))) 
      energy += P->TETRA_ENERGY[(ts - Tetraloops)/7];
  }
  return energy;
}

double interiorloop(int i, int j, int k, int l, int bp_type1, int bp_type2, int *intSequence) {
  double energy;
  int n1,n2;
  /* Interior loop, bulge or stack? */
  n1 = k-i-1;
  n2 = j-l-1;

  if ( (n1>0) && (n2>0) )
    /* Interior loop, special cases for interior loops of 
     * sizes 1+1, 1+2, 2+1, and 2+2. */
    if ( (n1==1) && (n2==1) )
      energy = P->int11[bp_type1][bp_type2]
	[intSequence[i+1]][intSequence[j-1]];
    else if ( (n1==2) && (n2==1) )
      energy = P->int21[bp_type2][bp_type1]
	[intSequence[j-1]][intSequence[i+1]][intSequence[i+2]];
    else if ( (n1==1) && (n2==2) )
      energy = P->int21[bp_type1][bp_type2]
	[intSequence[i+1]][intSequence[l+1]][intSequence[l+2]];
    else if ( (n1==2) && (n2==2) )
      energy = P->int22[bp_type1][bp_type2]
	  [intSequence[i+1]][intSequence[i+2]][intSequence[l+1]][intSequence[l+2]];
    else
      energy = ((n1+n2<=30)?(P->internal_loop[n1+n2]):
	(P->internal_loop[30]+(P->lxc*log((n1+n2)/30.))))+
	P->mismatchI[bp_type1][intSequence[i+1]][intSequence[j-1]] +
	P->mismatchI[bp_type2][intSequence[l+1]][intSequence[k-1]] +
	min2(MAX_NINIO, abs(n1-n2)*P->F_ninio[2]);
  else if ( (n1>0) || (n2>0) ) {
    /* Bulge */
    energy = ((n2+n1)<=30)?P->bulge[n1+n2]:
      (P->bulge[30]+(P->lxc*log((n1+n2)/30.)));
    if ( (n1+n2)==1 )
      /* A bulge of size one is so small that the base pairs
       * can stack */
      energy += P->stack[bp_type1][bp_type2];
    else {
      if ( bp_type1>2)
	energy += TerminalAU;
      if ( bp_type2>2)
	energy += TerminalAU;
    }
  }
  else { /* n1=n2=0 */
    /* Stacking base pair */
    energy = P->stack[bp_type1][bp_type2];
  }
  return energy;
}
