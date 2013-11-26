#ifndef DELTA_H
#define DELTA_H

#include <complex>

typedef std::complex<double> dcomplex;

void   neighbours(char *, int *);
int    numbp(int, int, int *);
int    bn(char);
void   initializeCanBasePair(int **);
void   translateToIntSequence(char *, int *);
void   initializeBasePairCounts(int **, int *, int );
double hairpinloop(int, int, int, short, short, char *);
double interiorloop(int, int, int, int, int, int, short, short, short, short);
void   solveSystem(dcomplex ***, char *, int *, int, int);
int    jPairedTo(int, int, int *);
int    jPairedIn(int, int, int *);
void   populateRemainingRoots(dcomplex ***, int, int, int);
void   populateMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, dcomplex ***, dcomplex *, int, int);
void   flushMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void   evaluateZ(int, dcomplex **, dcomplex **, dcomplex **, dcomplex **, dcomplex ***, dcomplex *, char *, char *, int *, int *, int **, int **, int, int, double);

#endif
