#ifndef DELTA_H
#define DELTA_H

#include "params.h"
#include <complex>

typedef std::complex<double> dcomplex;

void   neighbours(char *, int *);
void   pf(char *);
int    numbp(int, int, int *);
int    bn(char);
void   initializeCanBasePair(int[][5]);
void   translateToIntSequence(char *, int *);
void   initializeBasePairCounts(int **, int *, int );
double hairpinloop(int, int, int, int *, char *);
double interiorloop(int, int, int, int, int, int, int *);
double multiloop_closing(int, int, int, int, int, int, int *);
void   solveSystem(dcomplex ***, char *, int *, int, int);
int    jPairedTo(int, int, int *);
int    jPairedIn(int, int, int *);
void   populateRemainingRoots(dcomplex ***, int, int, int);
void   populateMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex ***, dcomplex *, int, int);
void   flushMatrices(dcomplex **, dcomplex **, dcomplex **, int);
void   evaluateZ(dcomplex **, dcomplex **, dcomplex **, dcomplex ***, dcomplex *, int, int);

#endif
