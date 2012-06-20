#ifndef DELTA_H
#define DELTA_H

#include "params.h"
#include <complex>

typedef std::complex<double> dcomplex;

/* To avoid computing MFE structures if this is not wanted. */
//#define COMPUTEMFE

/* Returns 1 if the two input characters are complementary bases */
int basepair(char , char );

void sub_structure(int , int , int *, int *);
void neighbours(char *,int *);
void pf(char *);
int basepaired_to(int ,int *);
int bp_diff(int *, int, int , int);
void print_bps(int *);

void solveSystem(int sequenceLength, dcomplex **rootsOfUnity, double *coefficients, double scalingFactor);
int numberOfBasePairs(int i, int j, int *basePairs);
int** fillBasePairCounts(int *basePairs, int n);
int jPairedTo(int i, int j, int *basePairs);
int jPairedIn(int i, int j, int *basePairs);
int canBasePair(int i, int j, char *sequence);
void printMatrix(dcomplex **matrix, char *title, int iStart, int iStop, int jStart, int jStop);

#endif
