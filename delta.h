#ifndef DELTA_H
#define DELTA_H

#include <complex>
#include <vector>
#include "memory_types.h"

// Use the dcomplex type from fftbor namespace (defined in memory_types.h)
// Provide a global alias for backward compatibility
using dcomplex = fftbor::dcomplex;

void neighbours(const char* inputSequence, const int* bpList);

int numbp(int i, int j, const int* bpList);
int bn(char A);

void initializeCanBasePair(fftbor::IntMatrix2D& canBasePair);
void translateToIntSequence(const char* a, int* intSequence);
void initializeBasePairCounts(fftbor::IntMatrix2D& numBasePairs, const int* bpList, int n);

double hairpinloop(int i, int j, int type, short si1, short sj1, const char* string);
double interiorloop(int i, int j, int k, int l, int type, int type_2, short si1, short sq1, short sj1, short sp1);

void solveSystem(fftbor::ComplexMatrix3D& solutions, const char* sequence, const int* structure, int sequenceLength, int runLength);

int jPairedTo(int i, int j, const int* basePairs);
int jPairedIn(int i, int j, const int* basePairs);

void populateRemainingRoots(fftbor::ComplexMatrix3D& solutions, int sequenceLength, int runLength, int lastRoot);

void populateMatrices(fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                      fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                      fftbor::ComplexMatrix3D& solutions, std::vector<dcomplex>& rootsOfUnity,
                      int sequenceLength, int runLength);

void flushMatrices(fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                   fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                   int sequenceLength);

void evaluateZ(int root,
               fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
               fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
               fftbor::ComplexMatrix3D& solutions, const std::vector<dcomplex>& rootsOfUnity,
               const char* inputSequence, const char* sequence,
               const int* intSequence, const int* bpList,
               const fftbor::IntMatrix2D& canBasePair, const fftbor::IntMatrix2D& numBasePairs,
               int sequenceLength, int runLength, double RT);

#endif
