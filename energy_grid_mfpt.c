#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "energy_grid_mfpt.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifdef __cplusplus
  extern "C" {
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);  
  }
#else
  extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
#endif

double** convertEnergyGridToTransitionMatrix(double* p, int length) {
  int i, j;
  double rowSum;
  double** transitionProbabilities = (double**)malloc(length * sizeof(double*));
  
  for (i = 0; i < length; ++i) {
    rowSum = 0.;
    
    transitionProbabilities[i] = (double*)malloc(length * sizeof(double));
      
    for (j = 0; j < length; ++j) {
      if (i != j) {
        transitionProbabilities[i][j] = MIN(
          1., 
          p[j] / p[i]
        ) / (length - 1);
        
        rowSum += transitionProbabilities[i][j];
      }
    }
    
    transitionProbabilities[i][i] = 1 - rowSum;
  }
  
  return transitionProbabilities;
}

double computeMFPT(int* k, int* l, double **transitionProbabilities, int length, int debug) {
  int i, j, x, y, startIndex, endIndex, inversionMatrixRowLength = length - 1;
  double mfptFromStart;
  
  for (i = 0, startIndex = -1, endIndex = -1; i < length; ++i) {
    if (k[i] == 0) {
      startIndex = i;
    }
    
    if (l[i] == 0) {
      endIndex = i;
    }
  }
  
  if (debug) {
    printf("startIndex:\t%d\n", startIndex);
    printf("endIndex:\t%d\n", endIndex);
  }
  
  if (startIndex < 0) {
    if (debug) {
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the starting state.\n");
    }
    return -1;
  }
  
  if (endIndex < 0) {
    if (debug) {
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the ending state.\n");
    }
    return -2;
  }
  
  double *mfpt            = (double*)calloc(inversionMatrixRowLength, sizeof(double));
  double *inversionMatrix = (double*)malloc((int)pow((double)inversionMatrixRowLength, 2.) * sizeof(double));
  
  // If startIndex > endIndex, we need to shift to the left by one because the endIndex row / column is being removed.
  if (startIndex > endIndex) {
    startIndex--;
  }
  
  if (debug) {
    printf("Inversion matrix:\n");
    printf("i\tj\tx\ty\tinversionMatrix[x, y]\n");
  }
  
  for (i = 0; i < length; ++i) {
    for (j = 0; j < length; ++j) { 
      if (i != endIndex && j != endIndex) {
        x = (i > endIndex ? i - 1 : i);
        y = (j > endIndex ? j - 1 : j);
        
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (endIndex) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversionMatrix and i, j are used for indexing into transitionProbabilities.
        inversionMatrix[x * inversionMatrixRowLength + y] = (i == j ? 1 - transitionProbabilities[i][j] : -transitionProbabilities[i][j]);
        
        if (debug) {
          printf("%d\t%d\t%d\t%d\t%f\n", i, j, x, y, inversionMatrix[x * inversionMatrixRowLength + y]);
        }
      }
    }
    
    if (debug) {
      printf("\n");
    }
  }
  
  inverse(inversionMatrix, inversionMatrixRowLength);
  
  for (i = 0; i < inversionMatrixRowLength; ++i) {
    for (j = 0; j < inversionMatrixRowLength; ++j) {
      mfpt[i] += inversionMatrix[i * inversionMatrixRowLength + j];
    }
  }
    
  mfptFromStart = mfpt[startIndex];
  free(mfpt);
  free(inversionMatrix);
  
  return mfptFromStart;
}

void inverse(double* A, int N) {
  // http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
  int *IPIV    = (int*)malloc((N + 1) * sizeof(int));
  int LWORK    = N * N;
  double* WORK = (double*)malloc(LWORK * sizeof(double));
  int INFO;

  dgetrf_(&N, &N, A, &N, IPIV, &INFO);
  dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

  free(IPIV);
  free(WORK);
}
