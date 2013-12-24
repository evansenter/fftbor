#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fftw3.h>
#include <limits>
#include "params.h"
#include "functions.h"
#include "shared/libmfpt_header.h"
#include "shared/libspectral_header.h"
#include "rna_misc_functions.h"
#include "vienna/functions.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

#define MIN_PAIR_DIST TURN
#define MAX_INTERIOR_DIST MAXLOOP
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

// #define SILENCE_OUTPUT  1
// #define TIMING_DEBUG    1
// #define FFTBOR_DEBUG    1
// #define TWIDDLE_DEBUG   1
// #define MEASURE_TWIDDLE 1
// #define OPENMP_DEBUG    1
// #define SINGLE_THREAD   1
// #define STRUCTURE_COUNT 1
#define DO_WORK

extern int    PRECISION, MAXTHREADS, ROW_LENGTH, MATRIX_FORMAT, SIMPLE_OUTPUT, TRANSITION_OUTPUT, SPECTRAL_OUTPUT, EXPLICIT_ENERGY_FILE, GLOBAL_SEQ_LENGTH;
extern double temperature;
extern char   *ENERGY, *GLOBAL_SEQ, *GLOBAL_STR_1, *GLOBAL_STR_2;
extern paramT *P;
double RT;
int    TWIDDLE, GLOBAL_BP_DIST;

extern "C" {
  void read_parameter_file(const char energyfile[]);
  int *get_iindx(unsigned int sequenceLength);
  unsigned int *maximumMatchingConstraint(const char *sequence, short *viennaBP);
}

void neighbors(char *inputSequence, int **bpList) {
  #ifdef TIMING_DEBUG
    struct timeval fullStart, fullStop, start, stop;
    gettimeofday(&fullStart, NULL);
    gettimeofday(&start, NULL);
  #endif
  
  int i, j, root, minimalRowLength, requestedRowLength, rowLength, runLength, numRoots, sequenceLength = strlen(inputSequence), inputStructureDist = 0, nonZeroCount = 0;
  double scalingFactor;
  RT      = 0.0019872370936902486 * (temperature + K0) * 100; // 0.01 * (kcal K) / mol
  TWIDDLE = 6;
  
  char precisionFormat[20];
  sprintf(precisionFormat, "%%+.0%df", PRECISION ? (int)floor(log(pow(2., PRECISION)) / log(10.)) : std::numeric_limits<double>::digits);

  char  *energyfile     = findEnergyFile();
  char  *sequence       = new char[sequenceLength + 1];
  short *intSequence    = (short *)xcalloc(sequenceLength + 1, sizeof(short));
  short **viennaBP      = new short*[2];
  int   ***numBasePairs = new int**[2];
  int   **canBasePair;
  int   *nonZeroIndices;
  
  // Load energy parameters.
  read_parameter_file(energyfile);
  P = scale_parameters();
  P -> model_details.special_hp = 1;

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
  GLOBAL_BP_DIST = inputStructureDist;
  
  // Secondary structure data structure in the slightly different format that Vienna uses.
  for (i = 0; i < 2; ++i) {
    viennaBP[i] = new short[sequenceLength + 1];
    
    viennaBP[i][0] = sequenceLength;
    for (j = 1; j <= sequenceLength; ++j) {
      viennaBP[i][j] = bpList[i][j] > 0 ? bpList[i][j] : 0;
    }
  }
  
  // Index for moving in quadratic distancy dimensions (from Vienna 2.1.2)
  int *index = get_iindx((unsigned)sequenceLength);
  
  // Maximally saturated structure constrained with the input structures.
  unsigned int *maxBPConstrained1 = maximumMatchingConstraint(inputSequence, viennaBP[0]);
  unsigned int *maxBPConstrained2 = maximumMatchingConstraint(inputSequence, viennaBP[1]);
  
  // Minimize the row size to the max BP distance.
  minimalRowLength = MAX2(
    bpList[0][0] + maxBPConstrained1[index[1] - sequenceLength], 
    bpList[1][0] + maxBPConstrained2[index[1] - sequenceLength]
  );
  
  // Initialize the row length and matrix size variables for the 2D evaluation.
  requestedRowLength = ROW_LENGTH > 0 ? (int)ceil(sequenceLength * ROW_LENGTH / 100.) : minimalRowLength;
  // Note: rowLength = (least even number >= sequenceLength) + 1 (for a seq. of length 9 rowLength = (0..10).length = 11)
  rowLength = (requestedRowLength % 2 ? requestedRowLength + 1 : requestedRowLength) + 1;
  // Note: runLength = (least number div. 4 >= rowLength ^ 2) / 2 (for a seq. of length 9 runLength = (11 ^ 2 + 3) / 2 = 62)
  runLength = ((int)pow((double)rowLength, 2) + ((int)pow((double)rowLength, 2) % 4)) / 2;
  numRoots  = runLength * 2;
  
  #if defined(_OPENMP) && defined(OPENMP_DEBUG)
    printf("Max threads possible: %d\n", omp_get_max_threads());
    
    #ifdef SINGLE_THREAD
      MAXTHREADS = 1;
    #endif
    
    printf("Setting number of threads: %d\n", MAXTHREADS);
  #endif
    
  #ifdef _OPENMP
    // Set number of threads for OpenMP
    omp_set_num_threads(MAXTHREADS);
  #endif
  
  // Initialize matricies for dynamic programming.
  dcomplex ***Z   = new dcomplex**[MAXTHREADS];
  dcomplex ***ZB  = new dcomplex**[MAXTHREADS];
  dcomplex ***ZM  = new dcomplex**[MAXTHREADS];
  dcomplex ***ZM1 = new dcomplex**[MAXTHREADS];

  #pragma omp parallel for default(shared)
  for (i = 0; i < MAXTHREADS; ++i) {
    Z[i]   = new dcomplex*[sequenceLength + 1];
    ZB[i]  = new dcomplex*[sequenceLength + 1];
    ZM[i]  = new dcomplex*[sequenceLength + 1];
    ZM1[i] = new dcomplex*[sequenceLength + 1];
    populateZMatrices(Z[i], ZB[i], ZM[i], ZM1[i], sequenceLength);
  }

  // Initialize convenience tables for storing solutions, roots of unity, probabilities and non-zero indices.
  dcomplex *solutions    = new dcomplex[runLength + 1];
  dcomplex *rootsOfUnity = new dcomplex[numRoots];
  double *probabilities  = (double *)xcalloc(2 * runLength + 1, sizeof(double));
  
  if (TRANSITION_OUTPUT || SPECTRAL_OUTPUT) {
    nonZeroIndices = (int *)xcalloc((int)pow(rowLength, 2.) + 1, sizeof(int));
  }
  
  // Populate convenience tables.
  populateMatrices(rootsOfUnity, numRoots);
  
  // Create convenience table for looking up powers of roots.
  dcomplex **rootToPower = new dcomplex*[MAXTHREADS];
  for (i = 0; i < MAXTHREADS; ++i) {
    rootToPower[i] = new dcomplex[(rowLength + 1) * sequenceLength + 1];
  }
  
  // Create convenience table for looking up 1D indexing of (k, l) coordinates.
  int **deltaTable = new int*[sequenceLength + 1];
  for (i = 0; i <= sequenceLength; ++i) {
    deltaTable[i] = new int[sequenceLength + 1];
    for (j = 0; j <= sequenceLength; ++j) {
      deltaTable[i][j] = DELTA_2D(i, j, rowLength);
    }
  } 
  
  // Create convenience table for boolean (i, j paired?) values.
  int **jPairedTo0 = new int*[sequenceLength + 1];
  int **jPairedTo1 = new int*[sequenceLength + 1];
  for (i = 0; i <= sequenceLength; ++i) {
    jPairedTo0[i] = new int[sequenceLength + 1];
    jPairedTo1[i] = new int[sequenceLength + 1];
    for (j = 0; j <= sequenceLength; ++j) {
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
    for (j = 0; j <= sequenceLength; ++j) {
      EM1[i][j] = new double[sequenceLength + 1];
      // Multiplied by a "twiddle" factor.
      EIL[i][j] = new double[TWIDDLE * (sequenceLength + 1)];
    }
  }
  
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "initialization")
    gettimeofday(&start, NULL);
  #endif
  
  calculateEnergies(
    inputSequence, 
    intSequence,
    canBasePair, 
    sequenceLength,
    RT, 
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
    printf("TWIDDLE:            %d\n", TWIDDLE);
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
      sequence,
      intSequence, 
      bpList, 
      canBasePair, 
      numBasePairs, 
      sequenceLength, 
      rowLength, 
      numRoots, 
      rootToPower[threadId], 
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

  solveSystem(probabilities, solutions, rowLength, runLength, inputStructureDist, nonZeroIndices, nonZeroCount, scalingFactor, precisionFormat);

  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "FFT")
    gettimeofday(&start, NULL);
  #endif
      
  #ifndef SILENCE_OUTPUT
    printOutput(probabilities, inputStructureDist, minimalRowLength, rowLength, nonZeroIndices, nonZeroCount, scalingFactor, precisionFormat);
  #endif
    
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "Print output / matrix inversion (if -X was provided)")
    gettimeofday(&start, NULL);
  #endif
    
  // Free memory.
  for (i = 0; i <= sequenceLength; ++i) {
    delete[] EH[i];
    delete[] EHM[i];
    delete[] EMA[i];
    delete[] EMB[i];
    delete[] EZ[i];
    for (j = 0; j <= sequenceLength; ++j) {
      delete[] EM1[i][j];
      delete[] EIL[i][j];
    }
    delete[] EM1[i];
    delete[] EIL[i];
  }
  delete[] EH;
  delete[] EHM;
  delete[] EMA;
  delete[] EMB;
  delete[] EZ;
  delete[] EIL;
  delete[] EM1;

  for (i = 0; i < 2; ++i) {
    for (j = 0; j < sequenceLength + 1; ++j) {
      free(numBasePairs[i][j]);
    }
    free(numBasePairs[i]);
  }

  delete[] numBasePairs;

  for (i = 0; i < 5; ++i) {
    free(canBasePair[i]);
  }
  free(canBasePair);
  free(intSequence);

  for (i = 0; i < MAXTHREADS; ++i) {
    for (j = 0; j <= sequenceLength; ++j) {
       delete[] Z[i][j];
       delete[] ZB[i][j];
       delete[] ZM[i][j];
       delete[] ZM1[i][j];
    }
    delete[] Z[i];
    delete[] ZB[i];
    delete[] ZM[i];
    delete[] ZM1[i];
  }

  delete[] Z;
  delete[] ZB;
  delete[] ZM;
  delete[] ZM1;

  for (i = 0; i < MAXTHREADS; ++i) {
    delete[] rootToPower[i];
  }
  delete[] rootToPower;

  for (i = 0; i <= sequenceLength; ++i) {
    delete[] deltaTable[i];
  }
  delete[] deltaTable;

  for (i = 0; i <= sequenceLength; ++i) {
    delete[] jPairedTo0[i];
    delete[] jPairedTo1[i];
  }

  delete[] jPairedTo0;
  delete[] jPairedTo1;
  
  delete[] nonZeroIndices;

  free(probabilities);
  delete[] solutions;
  delete[] rootsOfUnity;
  delete[] sequence; 
  
  #ifdef TIMING_DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "free memory")
  #endif
  
  #ifdef TIMING_DEBUG
    gettimeofday(&fullStop, NULL);
    TIMING(fullStart, fullStop, "total")
  #endif
}

void calculateEnergies(char *inputSequence, short *intSequence, int *canBasePair[5], int sequenceLength, double RT, double **EH, double ***EIL, double **EHM, double ***EM1, double **EMA, double **EMB, double **EZ) {
  int i, j, k, l, d, pos;
  
  #if MEASURE_TWIDDLE
    double maxTwiddle = 0;
  #endif

  for (i = 1; i <= sequenceLength; ++i) {
     for (d = MIN_PAIR_DIST + 1; d <= sequenceLength - i; ++d) {
      j = i + d;
      
      if (canBasePair[intSequence[i]][intSequence[j]]) {
        // ****************************************************************************
        // Solve ZB 
        // ****************************************************************************
        // In a hairpin, [i + 1, j - 1] unpaired.
        EH[i][j] = exp(-hairpinloop(i, j, canBasePair[intSequence[i]][intSequence[j]], intSequence[i + 1], intSequence[j - 1], inputSequence) / RT);

        pos = 0;
        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < MIN2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          // There can't be more than 30 unpaired bases between i and k, and there must be room between k and j for l
          for (l = MAX2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) { 
            // l needs to at least have room to pair with k, and there can be at most 30 unpaired bases between (i, k) + (l, j), with l < j
            if (canBasePair[intSequence[k]][intSequence[l]]) {
              #ifdef TWIDDLE_DEBUG
                if (pos >= TWIDDLE * (sequenceLength + 1)) {
                  fprintf(stderr, "Trying to access non-existant memory in EIL at index %d (%f).\n", pos, pos / (sequenceLength + 1.));
                }
              #endif
                
              #ifdef MEASURE_TWIDDLE
                maxTwiddle = pos / (sequenceLength + 1.) > maxTwiddle ? pos / (sequenceLength + 1.) : maxTwiddle;
              #endif
              
              // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
              EIL[i][j][pos] = exp(-interiorloop(i, j, k, l, canBasePair[intSequence[i]][intSequence[j]], canBasePair[intSequence[l]][intSequence[k]], intSequence[i + 1], intSequence[l + 1], intSequence[j - 1], intSequence[k - 1]) / RT);      
              pos++;
            }
          }
        }
        
        // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1)
        EHM[i][j] = exp(-(P->MLclosing + P->MLintern[canBasePair[intSequence[i]][intSequence[j]]]) / RT);
      }
      
      // ****************************************************************************
      // Solve ZM1
      // ****************************************************************************
      for (k = i + MIN_PAIR_DIST + 1; k <= j; ++k) {
        // k is the closing base pairing with i of a single component within the range [i, j]
        if (canBasePair[intSequence[i]][intSequence[k]]) {
          EM1[i][j][k] = exp(-(P->MLbase * (j - k) + P->MLintern[canBasePair[intSequence[i]][intSequence[k]]]) / RT);
        }
      }
        
      // ****************************************************************************
      // Solve ZM
      // ****************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        // Only one stem.
        EMA[i][k] = exp(-(P->MLbase * (k - i)) / RT);
        
        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
          EMB[j][k] = exp(-(P->MLintern[canBasePair[intSequence[k]][intSequence[j]]]) / RT);
        }
      }
        
      // **************************************************************************
      // Solve Z
      // **************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
        // (k, j) is the rightmost base pair in (i, j)
        if (canBasePair[intSequence[k]][intSequence[j]]) {
          EZ[j][k] = exp(-(canBasePair[intSequence[k]][intSequence[j]] > 2 ? P->TerminalAU : 0) / RT);

        }
      }
    }
  }
  
  #ifdef MEASURE_TWIDDLE
    printf("Max twiddle distance seen: %f\n", maxTwiddle);
  #endif
}

void evaluateZ(int root, dcomplex **Z, dcomplex **ZB, dcomplex **ZM, dcomplex **ZM1, dcomplex *solutions, dcomplex *rootsOfUnity, char *sequence, short *intSequence, int *bpList[2], int *canBasePair[5], int **numBasePairs[2], int sequenceLength, int rowLength, int numRoots, dcomplex *rootToPower, int **deltaTable, int **jPairedTo0, int **jPairedTo1, double **EH, double ***EIL, double **EHM, double ***EM1, double **EMA, double **EMB, double **EZ) {  
  int i, j, k, l, d, delta, pos;
  double energy;
  
  flushMatrices(Z, ZB, ZM, ZM1, sequenceLength);

  for (i = 0; i < (rowLength + 1) * sequenceLength + 1; ++i) {
    rootToPower[i] = ROOT_POW(root, i, numRoots);
  }
  
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
      
          ZB[i][j] += rootToPower[delta] * EH[i][j];
          
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
              
                ZB[i][j] += ZB[k][l] * rootToPower[delta] * EIL[i][j][pos];
                pos++;

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
         
            ZB[i][j] += ZM[i + 1][k - 1] * ZM1[k][j - 1] * rootToPower[delta] * energy;

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
          
            ZM1[i][j] += ZB[i][k] * EM1[i][j][k];
          
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
        
          ZM[i][j] += ZM1[k][j] * rootToPower[delta] * EMA[i][k];

          #ifdef STRUCTURE_COUNT
            ZM[j][i] += ZM1[j][k];
          #endif
        #endif

        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
          #ifdef DO_WORK
            delta  = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][i][k - 1] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][i][k - 1] - numBasePairs[1][k][j]];
          
            ZM[i][j] += ZM[i][k - 1] * ZM1[k][j] * rootToPower[delta] * EMB[j][k];

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
      
        Z[i][j] += Z[i][j - 1] * rootToPower[delta];

        #ifdef STRUCTURE_COUNT
          Z[j][i] += Z[j - 1][i];
        #endif
      #endif

      for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
        // (k, j) is the rightmost base pair in (i, j)
        if (canBasePair[intSequence[k]][intSequence[j]]) {
          #ifdef DO_WORK
            energy = EZ[j][k];
          #endif
            
          if (k == i) {
            #ifdef DO_WORK
              delta = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][k][j]];
            
              Z[i][j] += ZB[k][j] * rootToPower[delta] * energy;

              #ifdef STRUCTURE_COUNT
                Z[j][i] += ZB[j][k];
              #endif
            #endif
          } else {
            #ifdef DO_WORK
              delta = deltaTable[numBasePairs[0][i][j] - numBasePairs[0][i][k - 1] - numBasePairs[0][k][j]][numBasePairs[1][i][j] - numBasePairs[1][i][k - 1] - numBasePairs[1][k][j]];
            
              Z[i][j] += Z[i][k - 1] * ZB[k][j] * rootToPower[delta] * energy;

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
    printf(".");
  #endif
}

void solveSystem(double *probabilities, dcomplex *solutions, int rowLength, int runLength, int inputStructureDist, int *nonZeroIndices, int& nonZeroCount, double& scalingFactor, char *precisionFormat) {
  int i, x, y;
  double sum = 0, normalSum = 0;
  
  int solutionLength = (int)pow((double)rowLength, 2);
  int offset         = inputStructureDist % 2 ? 1 : 0;
  
  fftw_complex *signal = (fftw_complex *)malloc(runLength * sizeof(fftw_complex));
  fftw_complex *result = (fftw_complex *)malloc(runLength * sizeof(fftw_complex));
  
  fftw_plan plan = fftw_plan_dft_1d(runLength, signal, result, FFTW_BACKWARD, FFTW_ESTIMATE);
  sum            = 0;
  scalingFactor  = solutions[0].real();
  
  // For some reason it's much more numerically stable to set the signal via real / imag components separately.
  for (i = 0; i < runLength; ++i) {
    // Convert point-value solutions of Z(root) to 10^PRECISION * Z(root) / Z
    signal[i][FFTW_REAL] = pow(2., (double)PRECISION) * (solutions[i].real() / scalingFactor);
    signal[i][FFTW_IMAG] = pow(2., (double)PRECISION) * (solutions[i].imag() / scalingFactor);
  }
  
  #ifdef FFTBOR_DEBUG
    printf("Scaling factor: %f:\n", solutions[0].real());
    printf("Scaled solutions vector for the inverse DFT:\n");
    for (i = 0; i < runLength; ++i) {
      printf("%d: %+f %+fi\n", i, signal[i][FFTW_REAL], signal[i][FFTW_IMAG]);
    }
  #endif
  
  // Calculate transform, coefficients are in fftw_complex result array.
  fftw_execute(plan);
  
  for (i = 0; i < runLength; ++i) {
    // Truncate to user-specified precision; if set to 0, no truncation occurs (dangerous).
    if (!PRECISION) {
      solutions[i] = dcomplex(result[i][FFTW_REAL] / runLength, 0);
    } else {
      solutions[i] = dcomplex(pow(2., -PRECISION) * static_cast<int>(result[i][FFTW_REAL] / runLength), 0);
    }
    
    x = (2 * i + offset) / rowLength;
    y = (2 * i + offset) % rowLength;
    
    // Probabilities must be > 0 and satisfy the triangle inequality.
    if (
      solutions[i].real() > 0     && 
      x + y >= inputStructureDist &&
      x + inputStructureDist >= y &&
      y + inputStructureDist >= x
    ) {
      probabilities[2 * i + offset] = solutions[i].real();
      
      // Housekeeping.  
      sum += solutions[i].real();
      
      if (TRANSITION_OUTPUT || SPECTRAL_OUTPUT) {
        nonZeroIndices[nonZeroCount++] = 2 * i + offset;
      }
    }
  }
  
  // Normalizing pass.
  for (i = 0; i < solutionLength; ++i) {
    probabilities[i] /= sum;
    normalSum        += probabilities[i];
  }
  
  fftw_destroy_plan(plan);
  
  #ifdef FFTBOR_DEBUG
    printf("\nScaling factor: ");
    printf(precisionFormat, scalingFactor);
    printf("\nSum of eligible probabilities > 0: ");
    printf(precisionFormat, sum);
    printf("\nSum of normalized probabilities: ");
    printf(precisionFormat, normalSum);
    printf("\n");
  #endif
}

void printOutput(double *probabilities, int inputStructureDist, int minimalRowLength, int rowLength, int *nonZeroIndices, int& nonZeroCount, double& scalingFactor, char *precisionFormat) {
  int i, solutionLength = (int)pow((double)rowLength, 2);
  
  if (!(SIMPLE_OUTPUT || MATRIX_FORMAT || TRANSITION_OUTPUT || SPECTRAL_OUTPUT)) {
    printf("%d,%d,%d\nk\tl\tp(Z_{k,l}/Z)\t-RTln(Z_{k,l})\n", inputStructureDist, minimalRowLength, rowLength);
  }

  if (MATRIX_FORMAT) {
    for (i = 0; i < solutionLength; ++i) {
      if (i && !(i % rowLength)) {
        printf("\n");
      }
  
      printf(precisionFormat, probabilities[i]);
      printf("\t");
    }
    
    printf("\n");
  } else if (TRANSITION_OUTPUT) {
    calculateKinetics(nonZeroIndices, nonZeroCount, probabilities, rowLength, precisionFormat);
  } else if (SPECTRAL_OUTPUT) {
    populationProportion(nonZeroIndices, nonZeroCount, probabilities, rowLength, precisionFormat);
  } else {
    for (i = 0; i < solutionLength; ++i) { 
      if (probabilities[i] > 0) {
        printf("%d\t%d\t", i / rowLength, i % rowLength);
        printf(precisionFormat, probabilities[i]);
        printf("\t");
        printf(precisionFormat, -(RT / 100) * log(probabilities[i]) -(RT / 100) * log(scalingFactor));
        printf("\n");
      }
    }
  }
}

void calculateKinetics(int *nonZeroIndices, int nonZeroCount, double *probabilities, int rowLength, char *precisionFormat) {
  int error = 0;
  double mfpt;
  double** transitionMatrix;
  MFPT_PARAMETERS parameters;
  KLP_MATRIX klpMatrix;
  
  parameters                      = init_mfpt_params();  
  parameters.sequence_length      = GLOBAL_SEQ_LENGTH;
  parameters.distributed_epsilon  = 1e-8;
  parameters.bp_dist              = GLOBAL_BP_DIST;
  parameters.single_bp_moves_only = 1;
  parameters.hastings             = 1;  
  error                           = mfpt_error_handling(parameters);
  klpMatrix                       = init_klp_matrix(nonZeroCount);
  
  if (error) {
    fprintf(stderr, "Errors occured when calling the libmfpt files, terminating:\n");
    exit(0);
  }
  
  convert_fftbor2d_energy_grid_to_klp_matrix(nonZeroIndices, probabilities, rowLength, klpMatrix);
  transitionMatrix = convert_energy_grid_to_transition_matrix(&klpMatrix, parameters);
  mfpt             = compute_mfpt(klpMatrix, parameters, transitionMatrix);
  
  printf(precisionFormat, mfpt);
  printf("\n");
  
  free_transition_matrix(transitionMatrix, nonZeroCount);
  free_klp_matrix(klpMatrix);
}

void convert_fftbor2d_energy_grid_to_klp_matrix(int *nonZeroIndices, double *probabilities, int rowLength, KLP_MATRIX klpMatrix) {
  int i;
  
  for (i = 0; i < klpMatrix.length; ++i) {
    klpMatrix.k[i] = nonZeroIndices[i] / rowLength;
    klpMatrix.l[i] = nonZeroIndices[i] % rowLength;
    klpMatrix.p[i] = probabilities[nonZeroIndices[i]];
  }
}

void populationProportion(int *nonZeroIndices, int nonZeroCount, double *probabilities, int rowLength, char *precisionFormat) {
  int i, startIndex = -1, endIndex = -1, error = 0;
  double step_counter;
  double* transitionMatrix;
  EIGENSYSTEM eigensystem;
  SPECTRAL_PARAMS parameters;
  
  parameters = init_spectral_params();  
  error      = spectral_error_handling(parameters);
  
  if (error) {
    fprintf(stderr, "Errors occured when calling the libspectral files, terminating:\n");
    exit(0);
  }
  
  for (i = 0; i < nonZeroCount; ++i) {
    if (nonZeroIndices[i] / rowLength == 0 && nonZeroIndices[i] % rowLength == GLOBAL_BP_DIST) {
      startIndex = i;
    }
    
    if (nonZeroIndices[i] / rowLength == GLOBAL_BP_DIST && nonZeroIndices[i] % rowLength == 0) {
      endIndex = i;
    }
  }
  
  if (startIndex < 0 || endIndex < 0) {
    printf("The start index (%d) or the end index (%d) in the nonZeroIndices array doesn't exist, committing digital hari-kiri.\n", startIndex, endIndex);    
    printf("\n");
    printf("        /                    \n");
    printf("*//////{<>==================-\n");
    printf("        \\                   \n");
    exit(0);
  }
  
  transitionMatrix = convert_fftbor2d_energy_grid_to_transition_rate_matrix(nonZeroIndices, nonZeroCount, probabilities);
  eigensystem      = convert_transition_matrix_to_eigenvectors(transitionMatrix, nonZeroCount);
  invert_matrix(eigensystem);
  
  for (step_counter = parameters.start_time; step_counter <= parameters.end_time; step_counter += parameters.step_size) {
    printf(precisionFormat, step_counter);
    printf("\t");
    printf(precisionFormat, probability_at_time(eigensystem, pow(10, step_counter), startIndex, endIndex));
    printf("\t");
    printf(precisionFormat, probability_at_time(eigensystem, pow(10, step_counter), startIndex, startIndex));
    printf("\n");
  }

  free_eigensystem(eigensystem);
}

double* convert_fftbor2d_energy_grid_to_transition_rate_matrix(int *nonZeroIndices, int nonZeroCount, double *probabilities) {
  int i, j;
  double colSum;
  double* transitionMatrix;
  
  transitionMatrix = (double*)malloc(nonZeroCount * nonZeroCount * sizeof(double));
  
  for (i = 0; i < nonZeroCount; ++i) {
    colSum = 0;
    
    for (j = 0; j < nonZeroCount; ++j) {
      if (i != j) {
        transitionMatrix[i + nonZeroCount * j] = MIN2(1., probabilities[nonZeroIndices[j]] / probabilities[nonZeroIndices[i]]);
        colSum                                += transitionMatrix[i + nonZeroCount * j];
      }
      
      transitionMatrix[i + nonZeroCount * i] = -colSum;
    }
  }
  
  return transitionMatrix;
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

void populateMatrices(dcomplex *rootsOfUnity, int numRoots) {
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

// INLINE  PRIVATE int E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P){
inline double hairpinloop(int i, int j, int type, short si1, short sj1, char *string){
  double energy;
  int size = j-i-1;

  energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30]+(int)(P->lxc*log((size)/30.));
  if (P->model_details.special_hp){
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
        return (P->Tetraloop_E[(ts - P->Tetraloops)/7]);
    }
    else if (size == 6) {
      char tl[9]={0}, *ts;
      strncpy(tl, string, 8);
      if ((ts=strstr(P->Hexaloops, tl)))
        return (energy = P->Hexaloop_E[(ts - P->Hexaloops)/9]);
    }
    else if (size == 3) {
      char tl[6]={0,0,0,0,0,0}, *ts;
      strncpy(tl, string, 5);
      if ((ts=strstr(P->Triloops, tl))) {
        return (P->Triloop_E[(ts - P->Triloops)/6]);
      }
      return (energy + (type>2 ? P->TerminalAU : 0));
    }
  }
  energy += P->mismatchH[type][si1][sj1];

  return energy;
}

// INLINE  PRIVATE int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P){
inline double interiorloop(int i, int j, int k, int l, int type, int type_2, short si1, short sq1, short sj1, short sp1){
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns;
  double energy;
  int n1 = k-i-1;
  int n2 = j-l-1;
  energy = INF;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return P->stack[type][type_2];  /* stack */

  if (ns==0) {                      /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(P->internal_loop[nl+1]) : (P->internal_loop[30]+(int)(P->lxc*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1];}
      else if (nl==3){              /* 2x3 loop */
        energy = P->internal_loop[5]+P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]) : (P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

char* findEnergyFile() {
  if (EXPLICIT_ENERGY_FILE) {
    return ENERGY;
  } else {
    char* envPath, *tempPath, *splitPath, *energyLocation, *possiblePath;
    tempPath = getenv("PATH");
    envPath  = (char*)malloc((strlen(tempPath) + 2) * sizeof(char));
    strcpy(envPath, ".:");
    strcat(envPath, tempPath);
  
    if (envPath != NULL) {
      splitPath = strtok(envPath, ":");
    
      while (splitPath != NULL) {
        possiblePath = (char*)malloc((strlen(splitPath) + 1 + strlen(ENERGY)) * sizeof(char));
      
        strcpy(possiblePath, splitPath);
        strcat(possiblePath, "/");
        strcat(possiblePath, ENERGY);

        if (access(possiblePath, R_OK) != -1) {
          energyLocation = (char*)malloc(strlen(possiblePath) * sizeof(char));
          strcpy(energyLocation, possiblePath);
          splitPath = NULL;
        } else {
          splitPath = strtok(NULL, ":");
        }
      
        free(possiblePath);
      }
    }
    
    free(envPath);
    
    return (energyLocation != NULL ? energyLocation : ENERGY);
  }
}