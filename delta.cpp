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
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define FFTW_REAL 0
#define FFTW_IMAG 1
#define COMPLEX_CONJ(complexNumber) (dcomplex((complexNumber).real(), -(complexNumber).imag()))
#define DELTA_2D(expression1, expression2, n) ((expression1) * (n) + (expression2))
#define ROOT_POW(i, pow, n) (rootsOfUnity[((i) * (pow)) % (n)])
#define PRINT_COMPLEX(i, complex) printf("%d: %+f %+fi\n", i, complex[i].real(), complex[i].imag())
#define TIMING(start, stop, task) printf("Time in ms for %s: %.2f\n", task, (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0));
#define MIN2(A, B) ((A) < (B) ? (A) : (B))
#define MAX2(A, B) ((A) > (B) ? (A) : (B))

// #define SILENCE_OUTPUT 1
// #define TIMING_DEBUG 1
// #define FFTBOR_DEBUG 1
// #define OPENMP_DEBUG 1
// #define STRUCTURE_COUNT 1
// #define ENERGY_DEBUG (1 && !root)
#define DO_WORK

extern int    PRECISION, MAXTHREADS, ROW_LENGTH, MATRIX_FORMAT, SIMPLE_OUTPUT;
extern double temperature;
extern char   *ENERGY;
extern paramT *P;
double RT;

extern "C" void read_parameter_file(const char energyfile[]);

void neighbors(char *inputSequence, int **bpList) {
  #ifdef TIMING_DEBUG
    struct timeval fullStart, fullStop, start, stop;
    gettimeofday(&fullStart, NULL);
    gettimeofday(&start, NULL);
  #endif
  
  int i, root, requestedRowLength, rowLength, runLength, numRoots, sequenceLength = strlen(inputSequence), inputStructureDist = 0;
  RT = 0.0019872370936902486 * (temperature + 273.15) * 100; // 0.01 * (kcal K) / mol

  char *energyfile    = ENERGY;
  char *sequence      = new char[sequenceLength + 1];
  short *intSequence  = (short *)xcalloc(sequenceLength + 1, sizeof(short));
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
  for (i = 1; i <= sequenceLength; ++i) {
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
  
  // Initialize the row length and matrix size variables for the 2D evaluation.
  requestedRowLength = ROW_LENGTH > 0 ? ROW_LENGTH : sequenceLength;
  // Note: rowLength = (least even number >= sequenceLength) + 1 (for a seq. of length 9 rowLength = (0..10).length = 11)
  rowLength = (requestedRowLength % 2 ? requestedRowLength + 1 : requestedRowLength) + 1;
  // Note: runLength = (least number div. 4 >= rowLength ^ 2) / 2 (for a seq. of length 9 runLength = (11 ^ 2 + 3) / 2 = 62)
  runLength = ((int)pow((double)rowLength, 2) + ((int)pow((double)rowLength, 2) % 4)) / 2;
  numRoots  = runLength * 2;
  
  #if defined(_OPENMP) && defined(OPENMP_DEBUG)
    printf("Max threads possible: %d\n", omp_get_max_threads());
    printf("Setting number of threads: %d\n", MAXTHREADS);
  #endif
    
  // Set number of threads for OpenMP
  omp_set_num_threads(MAXTHREADS);
  
  // Initialize matricies for dynamic programming.
  dcomplex ***Z   = new dcomplex**[MAXTHREADS];
  dcomplex ***ZB  = new dcomplex**[MAXTHREADS];
  dcomplex ***ZM  = new dcomplex**[MAXTHREADS];
  dcomplex ***ZM1 = new dcomplex**[MAXTHREADS];

  #pragma omp parallel for default(shared)
  for (int threadId = 0; threadId < MAXTHREADS; ++threadId) {
    Z[threadId]   = new dcomplex*[sequenceLength + 1];
    ZB[threadId]  = new dcomplex*[sequenceLength + 1];
    ZM[threadId]  = new dcomplex*[sequenceLength + 1];
    ZM1[threadId] = new dcomplex*[sequenceLength + 1];
    populateZMatrices(Z[threadId], ZB[threadId], ZM[threadId], ZM1[threadId], sequenceLength);
  }

  // Initialize convenience tables for storing solutions, roots of unity.
  dcomplex *solutions    = new dcomplex[runLength + 1];
  dcomplex *rootsOfUnity = new dcomplex[numRoots];
  
  // Populate convenience tables.
  populateMatrices(solutions, rootsOfUnity, sequenceLength, numRoots);
  
  // Create convenience table for looking up powers of roots.
  dcomplex **rootToPower = new dcomplex*[runLength + 1];
  for (i = 0; i <= runLength; ++i) {
    rootToPower[i] = new dcomplex[(rowLength + 1) * sequenceLength + 1];
    for (int j = 0; j < (rowLength + 1) * sequenceLength + 1; ++j) {
      rootToPower[i][j] = ROOT_POW(i,j,numRoots);
    }
  }
  
  // Create convenience table for looking up 1D indexing of (k, l) coordinates.
  int **deltaTable = new int*[sequenceLength + 1];
  for (i = 0; i <= sequenceLength; ++i) {
    deltaTable[i] = new int[sequenceLength + 1];
    for (int j = 0; j <= sequenceLength; ++j) {
      deltaTable[i][j] = DELTA_2D(i,j,rowLength);
    }
  } 
  
  // Create convenience table for boolean (i, j paired?) values.
  int **jPairedTo0 = new int*[sequenceLength + 1];
  int **jPairedTo1 = new int*[sequenceLength + 1];
  for (i = 0; i <= sequenceLength; ++i) {
    jPairedTo0[i] = new int[sequenceLength + 1];
    jPairedTo1[i] = new int[sequenceLength + 1];
    for (int j = 0; j <= sequenceLength; ++j) {
      jPairedTo0[i][j] = jPairedTo(i,j,bpList[0]);
      jPairedTo1[i][j] = jPairedTo(i,j,bpList[1]);
    }
  }

  // Initialize tables for precalculating energies.
  double **EZ   = new double*[sequenceLength + 1];
  double **EH   = new double*[sequenceLength + 1];
  double **EHM  = new double*[sequenceLength + 1];
  double **EMA  = new double*[sequenceLength + 1]; 
  double **EMB  = new double*[sequenceLength + 1];
  double ***EIL = new double**[sequenceLength + 1];
  double ***EM1 = new double**[sequenceLength + 1]; 
  
  
  for (i = 0; i <= sequenceLength; ++i) {
    EZ[i]  = new double[sequenceLength + 1];
    EH[i]  = new double[sequenceLength + 1];
    EHM[i] = new double[sequenceLength + 1];
    EMA[i] = new double[sequenceLength + 1];
    EMB[i] = new double[sequenceLength + 1];
    EIL[i] = new double*[sequenceLength + 1];
    EM1[i] = new double*[sequenceLength + 1];
    for (int j = 0; j <= sequenceLength; ++j) {
      EM1[i][j] = new double[sequenceLength + 1];
      // 5 is still a "twiddle" factor.
      EIL[i][j] = new double[5 * (sequenceLength + 1)];
    }
  }
  
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "initialization")
    gettimeofday(&start, NULL);
  #endif
  
  calculateEnergies(
    inputSequence, 
    sequence, 
    intSequence,
    bpList,
    canBasePair, 
    numBasePairs, 
    sequenceLength,
    RT, 
    rootToPower,
    deltaTable,
    jPairedTo0,
    jPairedTo1,
    EH,
    EIL,
    EHM,
    EM1,
    EMA,
    EMB,
    EZ
  );
  
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "pre-calculate energies")
    gettimeofday(&start, NULL);
  #endif
  
  #ifdef FFTBOR_DEBUG
    printf("sequenceLength:     %d\n", sequenceLength);
    printf("rowLength:          %d\n", rowLength);
    printf("runLength:          %d\n", runLength);
    printf("numRoots:           %d\n", numRoots);
    printf("inputStructureDist: %d\n", inputStructureDist);
    printf("Roots of unity:\n");
    for (root = 0; root < numRoots; ++root) {
      PRINT_COMPLEX(root, rootsOfUnity);
    }
  #endif

  // Start main recursions (root <= maxroot is an optimization for roots of unity).
  int maxroot = runLength / 2;
  int threadId;
    
  #pragma omp parallel for private(root, threadId) shared(Z, ZB, ZM, ZM1, bpList, canBasePair, numBasePairs, inputStructureDist, sequenceLength, runLength, rowLength, numRoots, RT, maxroot, intSequence, sequence, inputSequence, rootsOfUnity, solutions, rootToPower, deltaTable, jPairedTo0, jPairedTo1, EH, EIL, EHM, EM1, EMA, EMB, EZ) default(none) num_threads(MAXTHREADS)
  for (root = 0; root <= maxroot; ++root) {
    #ifdef _OPENMP
      threadId = omp_get_thread_num();
    #else
      threadId = 0;
    #endif
    
    evaluateZ(
      root, 
      Z[threadId], 
      ZB[threadId], 
      ZM[threadId], 
      ZM1[threadId], 
      solutions, 
      rootsOfUnity, 
      inputSequence, 
      sequence, 
      intSequence, 
      bpList, 
      canBasePair, 
      numBasePairs, 
      inputStructureDist, 
      sequenceLength, 
      rowLength, 
      numRoots, 
      RT, 
      rootToPower, 
      deltaTable, 
      jPairedTo0, 
      jPairedTo1,
      EH,
      EIL,
      EHM,
      EM1,
      EMA,
      EMB,
      EZ
    );
  }
  
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "explicit evaluation of Z")
    gettimeofday(&start, NULL);
  #endif

  // Convert point-value solutions to coefficient form w/ inverse DFT.
  populateRemainingRoots(solutions, rootsOfUnity, runLength, inputStructureDist);
  
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "inferred evaluation of Z")
    gettimeofday(&start, NULL);
  #endif
  
  // Free memory.
  for (i = 0; i <= sequenceLength; ++i) {
    delete [] EH[i];
    delete [] EHM[i];
    delete [] EMA[i];
    delete [] EMB[i];
    delete [] EZ[i];
    for (int j = 0; j <= sequenceLength; ++j) {
      delete [] EM1[i][j];
      delete [] EIL[i][j];
    }
    delete [] EM1[i];
    delete [] EIL[i];
  }
  delete [] EH;
  delete [] EHM;
  delete [] EMA;
  delete [] EMB;
  delete [] EZ;
  delete [] EIL;
  delete [] EM1;
  for (i = 0; i < 2; ++i) {
    for (int j = 0; j < sequenceLength + 1; ++j) {
      free(numBasePairs[i][j]);
    }
    free(numBasePairs[i]);
  }
  
  delete [] numBasePairs;

  for (i = 0; i < 5; ++i) {
    free(canBasePair[i]);
  }
  free(canBasePair);
  free(intSequence);

  for (int threadId = 0; threadId < MAXTHREADS; ++threadId) {
    for (int j = 0; j <= sequenceLength; ++j) {
       delete [] Z[threadId][j];
       delete [] ZB[threadId][j];
       delete [] ZM[threadId][j];
       delete [] ZM1[threadId][j];
    }
    delete [] Z[threadId];
    delete [] ZB[threadId];
    delete [] ZM[threadId];
    delete [] ZM1[threadId];
  }

  delete [] Z;
  delete [] ZB;
  delete [] ZM;
  delete [] ZM1;

  for (i = 0; i <= runLength; ++i) {
    delete [] rootToPower[i];
  }
  delete [] rootToPower;

  for (i = 0; i <= sequenceLength; ++i) {
    delete [] deltaTable[i];
  }
  delete [] deltaTable;

  for (i = 0; i <= sequenceLength; ++i) {
    delete [] jPairedTo0[i];
    delete [] jPairedTo1[i];
  }

  delete [] jPairedTo0;
  delete [] jPairedTo1;
 
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "free memory")
    gettimeofday(&start, NULL);
  #endif

  solveSystem(solutions, rootsOfUnity, sequence, bpList, sequenceLength, rowLength, runLength, inputStructureDist);

  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "FFT")
  #endif

  delete [] solutions;
  delete [] rootsOfUnity;
  delete [] sequence; 
  
  #ifdef TIMING_DEBUG
    gettimeofday(&fullStop, NULL);
    TIMING(fullStart, fullStop, "total")
  #endif
}

void calculateEnergies(char *inputSequence, char *sequence, short *intSequence, int *bpList[2], int *canBasePair[5], int **numBasePairs[2], int sequenceLength, double RT, dcomplex **rootToPower, int **deltaTable, int **jPairedTo0, int **jPairedTo1, double **EH, double ***EIL, double **EHM, double ***EM1, double **EMA, double **EMB, double **EZ) {
  int i, j, k, l, d, delta, pos;

  #pragma opm parallel for shared(default)  
  for (i = 1; i <= sequenceLength; ++i) {
     
     for (d = MIN_PAIR_DIST + 1; d <= sequenceLength-i; ++d) {
      j = i + d;
      
      if (canBasePair[intSequence[i]][intSequence[j]]) {
        // ****************************************************************************
        // Solve ZB 
        // ****************************************************************************
          // In a hairpin, [i + 1, j - 1] unpaired.
           EH[i][j]  = exp(-hairpinloop(i, j, canBasePair[intSequence[i]][intSequence[j]], intSequence[i+1], intSequence[j-1], inputSequence)/RT);

        pos=0;
        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < MIN2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          // There can't be more than 30 unpaired bases between i and k, and there must be room between k and j for l
          for (l = MAX2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) { 
            // l needs to at least have room to pair with k, and there can be at most 30 unpaired bases between (i, k) + (l, j), with l < j
            if (canBasePair[intSequence[k]][intSequence[l]]) {
                // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
                EIL[i][j][pos] = exp(-interiorloop(i, j, k, l, canBasePair[intSequence[i]][intSequence[j]], canBasePair[intSequence[l]][intSequence[k]], intSequence[i+1], intSequence[l+1], intSequence[j-1], intSequence[k-1])/RT);      
                pos++;
            }
          }
        }
        
        // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1)
        EHM[i][j] = exp(-(P->MLclosing + P->MLintern[canBasePair[intSequence[i]][intSequence[j]]])/RT);
      }
      
      // ****************************************************************************
      // Solve ZM1
      // ****************************************************************************
      for (k = i + MIN_PAIR_DIST + 1; k <= j; ++k) {
        // k is the closing base pairing with i of a single component within the range [i, j]
        if (canBasePair[intSequence[i]][intSequence[k]]) {
            EM1[i][j][k] = exp(-(P->MLbase * (j - k) + P->MLintern[canBasePair[intSequence[i]][intSequence[k]]])/RT);
        }
      }
        
      // ****************************************************************************
      // Solve ZM
      // ****************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
          // Only one stem.
          EMA[i][k] = exp(-(P->MLbase * (k - i))/RT);
        
        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
            EMB[j][k] = exp(-(P->MLintern[canBasePair[intSequence[k]][intSequence[j]]])/RT);
        }
      }
        
      // **************************************************************************
      // Solve Z
      // **************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
        // (k, j) is the rightmost base pair in (i, j)
        if (canBasePair[intSequence[k]][intSequence[j]]) {
            EZ[j][k] = exp(-(canBasePair[intSequence[k]][intSequence[j]] > 2 ? TerminalAU37 : 0)/RT);

        }
      }
    }
  }
}

void evaluateZ(int root, dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, dcomplex *solutions, dcomplex *rootsOfUnity, char *inputSequence, char *sequence, short *intSequence, int *bpList[2], int *canBasePair[5], int **numBasePairs[2], int inputStructureDist, int sequenceLength, int rowLength, int numRoots, double RT, dcomplex **rootToPower, int **deltaTable, int **jPairedTo0, int **jPairedTo1,double **EH, double ***EIL, double **EHM, double ***EM1, double **EMA, double **EMB, double **EZ) {  
  int i, j, k, l, d, delta, pos;
  double energy;
  
  flushMatrices(Z, ZB, ZM, ZM1, sequenceLength);
    
  #ifdef ENERGY_DEBUG
    printf("RT: %f\n", RT / 100);
  #endif
  
  for (d = MIN_PAIR_DIST + 1; d < sequenceLength; ++d) {
    for (i = 1; i <= sequenceLength - d; ++i) {
      j = i + d;
      
      if (canBasePair[intSequence[i]][intSequence[j]]) {
        // ****************************************************************************
        // Solve ZB 
        // ****************************************************************************
        #ifdef DO_WORK
          // In a hairpin, [i + 1, j - 1] unpaired.
          delta  = deltaTable[numBasePairs[0][i][j] + jPairedTo0[i][j]][numBasePairs[1][i][j] + jPairedTo1[i][j]];
          
      
          ZB[i][j] += rootToPower[root][delta] * EH[i][j];//exp(-energy / RT);

          #ifdef ENERGY_DEBUG
            printf("%+f: GetHairpinEnergy(%c (%d), %c (%d)); Delta = %d\n", energy / 100, sequence[i], i, sequence[j], j, delta);
          #endif

          #ifdef STRUCTURE_COUNT
            ZB[j][i] += 1;
          #endif
        #endif
        
        pos = 0;
        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < MIN2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          // There can't be more than 30 unpaired bases between i and k, and there must be room between k and j for l
          for (l = MAX2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) { 
            // l needs to at least have room to pair with k, and there can be at most 30 unpaired bases between (i, k) + (l, j), with l < j
            if (canBasePair[intSequence[k]][intSequence[l]]) {
              #ifdef DO_WORK
                // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
                delta  = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][k][l] + jPairedTo0[i][j]][numBasePairs[1][i][j] - numBasePairs[1][k][l] + jPairedTo1[i][j]];
              
                ZB[i][j] += ZB[k][l] * rootToPower[root][delta] * EIL[i][j][pos];//exp(-energy / RT);
                pos++;
                #ifdef ENERGY_DEBUG
                  printf("%+f: GetInteriorStackingAndBulgeEnergy(%c (%d), %c (%d), %c (%d), %c (%d)); Delta = %d\n", energy / 100, sequence[i], i, sequence[j], j, sequence[k], k, sequence[l], l, delta);
                #endif

                #ifdef STRUCTURE_COUNT
                  ZB[j][i] += ZB[l][k];
                #endif
              #endif
            }
          }
        }
        energy = EHM[i][j];
        for (k = i + MIN_PAIR_DIST + 3; k < j - MIN_PAIR_DIST - 1; ++k) {
          #ifdef DO_WORK
            // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1)
            delta  = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][i + 1][k - 1] - numBasePairs[0][k][j - 1] + jPairedTo0[i][j]][numBasePairs[1][i][j] - numBasePairs[1][i + 1][k - 1] - numBasePairs[1][k][j - 1] + jPairedTo1[i][j]];
         
            ZB[i][j] += ZM[i + 1][k - 1] * ZM1[k][j - 1] * rootToPower[root][delta] * energy;//exp(-energy / RT);
                 
            #ifdef ENERGY_DEBUG
              printf("%+f: MultiloopA + MultiloopB; Delta = %d\n", energy / 100, delta);
            #endif

            #ifdef STRUCTURE_COUNT
              ZB[j][i] += ZM[k - 1][i + 1] * ZM1[j - 1][k];
            #endif
          #endif
        }
      }
      
      // ****************************************************************************
      // Solve ZM1
      // ****************************************************************************
      for (k = i + MIN_PAIR_DIST + 1; k <= j; ++k) {
        // k is the closing base pairing with i of a single component within the range [i, j]
        if (canBasePair[intSequence[i]][intSequence[k]]) {
          #ifdef DO_WORK
            delta  = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][i][k]][numBasePairs[1][i][j] - numBasePairs[1][i][k]];
          
            ZM1[i][j] += ZB[i][k] * EM1[i][j][k];//exp(-energy / RT);
          
            #ifdef ENERGY_DEBUG
              printf("%+f: MultiloopB + MultiloopC * (%d - %d); Delta = %d\n", energy / 100, j, k, delta);
            #endif
          
            #ifdef STRUCTURE_COUNT
              ZM1[j][i] += ZB[k][i];
            #endif
          #endif
        }
      }
        
      // ****************************************************************************
      // Solve ZM
      // ****************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        #ifdef DO_WORK
          // Only one stem.
          delta  = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][k][j]];
        
          ZM[i][j] += ZM1[k][j] * rootToPower[root][delta] * EMA[i][k];//exp(-energy / RT);

          #ifdef ENERGY_DEBUG
            printf("%+f: MultiloopC * (%d - %d); Delta = %d\n", energy / 100, k, i, delta);
          #endif

          #ifdef STRUCTURE_COUNT
            ZM[j][i] += ZM1[j][k];
          #endif
        #endif

        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
          #ifdef DO_WORK
            delta  = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][i][k - 1] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][i][k - 1] - numBasePairs[1][k][j]];
          
            ZM[i][j] += ZM[i][k - 1] * ZM1[k][j] * rootToPower[root][delta] * EMB[j][k];//exp(-energy / RT);

            #ifdef ENERGY_DEBUG
              printf("%+f: MultiloopB; Delta = %d\n", energy / 100, delta);
            #endif

            #ifdef STRUCTURE_COUNT
              ZM[j][i] += ZM[k - 1][i] * ZM1[j][k];
            #endif
          #endif
        }
      }
        
      // **************************************************************************
      // Solve Z
      // **************************************************************************
      #ifdef DO_WORK
        delta = deltaTable[jPairedIn(i, j, bpList[0])][jPairedIn(i, j, bpList[1])];
      
        Z[i][j] += Z[i][j - 1] * rootToPower[root][delta];

        #ifdef STRUCTURE_COUNT
          Z[j][i] += Z[j - 1][i];
        #endif
      #endif

      for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
        // (k, j) is the rightmost base pair in (i, j)
        if (canBasePair[intSequence[k]][intSequence[j]]) {
          #ifdef DO_WORK
            energy = EZ[j][k];//canBasePair[intSequence[k]][intSequence[j]] > 2 ? TerminalAU37 : 0;

            #ifdef ENERGY_DEBUG
              printf("%+f: %c-%c == (2 || 3) ? 0 : GUAU_penalty; Delta = %d\n", energy / 100, sequence[k], sequence[j], delta);
            #endif
          #endif
            
          if (k == i) {
            #ifdef DO_WORK
              delta = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][k][j]];
            
              Z[i][j] += ZB[k][j] * rootToPower[root][delta] * energy;//exp(-energy / RT);

              #ifdef STRUCTURE_COUNT
                Z[j][i] += ZB[j][k];
              #endif
            #endif
          } else {
            #ifdef DO_WORK
              delta = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][i][k - 1] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][i][k - 1] - numBasePairs[1][k][j]];
            
              Z[i][j] += Z[i][k - 1] * ZB[k][j] * rootToPower[root][delta] * energy;//exp(-energy / RT);

              #ifdef STRUCTURE_COUNT
                Z[j][i] += Z[k - 1][i] * ZB[j][k];
              #endif
            #endif
          }
        }
      }
    }
  }
  
  solutions[root] = Z[1][sequenceLength];

  #ifdef FFTBOR_DEBUG
    std::cout << '.' << std::flush;
  #endif
}

void solveSystem(dcomplex *solutions, dcomplex *rootsOfUnity, char *sequence, int **structure, int sequenceLength, int rowLength, int runLength, int inputStructureDist) {
  char precisionFormat[20];
  int i, solutionLength = (int)pow((double)rowLength, 2);
  double *probabilities = (double *)xcalloc(2*runLength+1, sizeof(double));
  double scalingFactor, sum;
  
  sprintf(precisionFormat, "%%+.0%df", PRECISION ? PRECISION : std::numeric_limits<double>::digits10);
  
  fftw_complex signal[runLength];
  fftw_complex result[runLength];
  
  fftw_plan plan = fftw_plan_dft_1d(runLength, signal, result, FFTW_BACKWARD, FFTW_ESTIMATE);
  sum            = 0;
  scalingFactor  = solutions[0].real();
  
  // For some reason it's much more numerically stable to set the signal via real / imag components separately.
  for (i = 0; i < runLength; ++i) {
    // Convert point-value solutions of Z(root) to 10^PRECISION * Z(root) / Z
    signal[i][FFTW_REAL] = (pow(10.0, (double)PRECISION) * solutions[i].real()) / scalingFactor;
    signal[i][FFTW_IMAG] = (pow(10.0, (double)PRECISION) * solutions[i].imag()) / scalingFactor;
  }
  
  #ifdef FFTBOR_DEBUG    
    std::cout << "Scaling factor: " << solutions[0] << std::endl;
    printf("Scaled solutions vector for the inverse DFT:\n");
    for (i = 0; i < runLength; ++i) {
      printf("%d: %+f %+fi\n", i, signal[i][FFTW_REAL], signal[i][FFTW_IMAG]);
    }
  #endif
  
  // Calculate transform, coefficients are in fftw_complex result array.
  fftw_execute(plan);
  
  for (i = 0; i < runLength; ++i) {
    // Truncate to user-specified precision, default is 4 and if set to 0, no truncation occurs (dangerous).
    if (!PRECISION) {
      solutions[i] = dcomplex(result[i][FFTW_REAL] / runLength, 0);
    } else {
      solutions[i] = dcomplex(pow(10.0, -PRECISION) * static_cast<int>(result[i][FFTW_REAL] / runLength), 0);
    }
    
    probabilities[inputStructureDist % 2 ? 2 * i + 1 : 2 * i] = solutions[i].real();
        
    sum += solutions[i].real();
  }
  
  #ifndef SILENCE_OUTPUT
    if (SIMPLE_OUTPUT) {
      for (i = 0; i < solutionLength; ++i) { 
        if (probabilities[i] > 0) {
          printf("%d\t%d\t", i / rowLength, i % rowLength);
          printf(precisionFormat, probabilities[i]);
          printf("\t%f\n", -(RT / 100) * log(probabilities[i]) -(RT / 100) * log(scalingFactor));
        }
      }
    } else if (MATRIX_FORMAT) {
      for (i = 0; i < solutionLength; ++i) {
        if (!(i % rowLength)) {
          printf("\n");
        }
    
        printf(precisionFormat, probabilities[i]);
        printf("\t");
      }
      
      printf("\n");
    } else {
      printf("\nk\tl\tp(Z_{k,l}/Z)\t-RTln(Z_{k,l})\n");
      
      for (i = 0; i < solutionLength; ++i) { 
        printf("%d\t%d\t", i / rowLength, i % rowLength);
        printf(precisionFormat, probabilities[i]);
        printf("\t%f\n", -(RT / 100) * log(probabilities[i]) -(RT / 100) * log(scalingFactor));
      }
    }
    
    #ifdef FFTBOR_DEBUG
      printf("\nScaling factor: %.15f\n", scalingFactor);
      std::cout << "Sum: " << sum << std::endl << std::endl;
    #endif
  #endif
  
  fftw_destroy_plan(plan);
  free(probabilities);
}

inline int jPairedTo(int i, int j, int *basePairs) {
  return basePairs[i] == j ? -1 : 1;
}

inline int jPairedIn(int i, int j, int *basePairs) {
  return basePairs[j] >= i && basePairs[j] < j ? 1 : 0;
}

void populateRemainingRoots(dcomplex *solutions, dcomplex *rootsOfUnity, int runLength, int inputStructureDist) {
  // Optimization leveraging complex conjugate of solutions to the polynomial.
  int i;
  
  #ifdef FFTBOR_DEBUG    
    printf("Solutions (before populating remaining roots):\n");
    for (i = 0; i < runLength; ++i) {
      PRINT_COMPLEX(i, solutions);
    }
  #endif
 
  for (i = runLength / 2 + 1; i < runLength; ++i) {
    #ifdef FFTBOR_DEBUG
      printf("l: %d, r %d\n", i, runLength - i);
    #endif
    solutions[i] = COMPLEX_CONJ(solutions[runLength - i]);
  }
  
  #ifdef FFTBOR_DEBUG    
    printf("Solutions (after populating remaining roots):\n");
    for (i = 0; i < runLength; ++i) {
      PRINT_COMPLEX(i, solutions);
    }
  #endif
  
  if (inputStructureDist % 2) {
    for (i = 0; i < runLength; ++i) {
      solutions[i] = COMPLEX_CONJ(rootsOfUnity[i]) * solutions[i];
    }
    for (i = runLength / 2 + 1; i <= runLength; ++i) {
      solutions[i] = dcomplex(-1, 0) * solutions[i];
    }
    
    #ifdef FFTBOR_DEBUG
      printf("Solutions (after multiplying solutions by nu^{-k} (and the scalar -1 for solutions y_{k: k > M_{0} / 2})):\n");
      for (i = 0; i < runLength; ++i) {
        PRINT_COMPLEX(i, solutions);
      }
    #endif
  }
}

void populateZMatrices(dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, int sequenceLength){
  int i;
  for (i = 0; i <= sequenceLength; ++i) {
    Z[i]   = new dcomplex[sequenceLength + 1];
    ZB[i]  = new dcomplex[sequenceLength + 1];
    ZM[i]  = new dcomplex[sequenceLength + 1];
    ZM1[i] = new dcomplex[sequenceLength + 1];
  }
}

void populateMatrices(dcomplex *solutions, dcomplex *rootsOfUnity, int sequenceLength, int numRoots) {
  int i;
  #pragma omp parallel for default(shared)
  for (i = 0; i < numRoots; ++i) {
    rootsOfUnity[i] = dcomplex(cos(-2 * M_PI * i / numRoots), sin(-2 * M_PI * i / numRoots));
  }
}

inline void flushMatrices(dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, int sequenceLength) {
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

short bn(char A) {
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

void translateToIntSequence(char *a, short *intSequence) {
  int i;
  intSequence[0] = strlen(a);
  for (i=0; i < intSequence[0]; ++i)
    intSequence[i+1] = bn(a[i]);
}

void initializeBasePairCounts(int **numBasePairs, int *bpList, int n) {
  int d, i, j;
  for (i = 1; i <= n; ++i)
    for (j = 1; j <= n; ++j)
      numBasePairs[i][j] = 0;
  for (d = MIN_PAIR_DIST+1; d < n; d++) {
    for (i = 1; i <= n - d; ++i) {
      j = i + d;
      numBasePairs[i][j] = numbp(i, j, bpList);
    }
  }
}

inline double hairpinloop(int i, int j, int bp_type, short iplus1, short jminus1, char *a) {
  // char *a is intended to be 0-indexed.
  double energy;
  energy = ((j-i-1) <= 30) ? P->hairpin[j-i-1] :
    P->hairpin[30]+(P->lxc*log((j-i-1)/30.));
  if ((j-i-1) >3) /* No mismatch for triloops */
    energy += P->mismatchH[bp_type][iplus1][jminus1];
  else {
    char tl[6]={0}, *ts;
    strncpy(tl, a+i-1, 5);
    if ((ts=strstr(Triloops, tl))) 
      energy += P->Triloop_E[(ts - Triloops)/6];
    if (bp_type>2)
      energy += TerminalAU37;
	}
  if ((j-i-1) == 4) { /* check for tetraloop bonus */
    char tl[7]={0}, *ts;
    strncpy(tl, a+i-1, 6);
    if ((ts=strstr(Tetraloops, tl))) 
      energy += P->TETRA_ENERGY[(ts - Tetraloops)/7];
  }
  return energy;
}

inline double interiorloop(int i, int j, int k, int l, int bp_type1, int bp_type2, short iplus1, short lplus1, short jminus1, short kminus1) {
  double energy;
  int n1,n2;
  /* Interior loop, bulge or stack? */
  n1 = k-i-1;
  n2 = j-l-1;

  if ( (n1>0) && (n2>0) ){
    if(n1>2 || n2>2)
      return ((n1+n2<=30)?(P->internal_loop[n1+n2]):
        (P->internal_loop[30]+(P->lxc*log((n1+n2)/30.))))+
        P->mismatchI[bp_type1][iplus1][jminus1] +
        P->mismatchI[bp_type2][lplus1][kminus1] +
        MIN2(MAX_NINIO, abs(n1-n2)*P->F_ninio[2]);
    /* Interior loop, special cases for interior loops of 
     * sizes 1+1, 1+2, 2+1, and 2+2. */
    
    else if ( (n1==1) && (n2==1) )
      return P->int11[bp_type1][bp_type2][iplus1][jminus1];
    else if ( (n1==2) && (n2==1) )
      return P->int21[bp_type2][bp_type1][lplus1][iplus1][kminus1];
    else if ( (n1==1) && (n2==2) )
      return P->int21[bp_type1][bp_type2][iplus1][lplus1][jminus1];
    else //if ( (n1==2) && (n2==2) )
      return P->int22[bp_type1][bp_type2][iplus1][kminus1][lplus1][jminus1];
    
  }else if ( (n1>0) || (n2>0) ) {
    /* Bulge */
    energy = ((n2+n1)<=30)?P->bulge[n1+n2]:
      (P->bulge[30]+(P->lxc*log((n1+n2)/30.)));
    if ( (n1+n2)==1 )
      /* A bulge of size one is so small that the base pairs
       * can stack */
      return energy + P->stack[bp_type1][bp_type2];
    else {
      if ( bp_type1>2)
	energy += TerminalAU37;
      if ( bp_type2>2)
	energy += TerminalAU37;
      return energy;
    }
  }
  else { /* n1=n2=0 */
    /* Stacking base pair */
    return P->stack[bp_type1][bp_type2];
  }
  return energy;
}


