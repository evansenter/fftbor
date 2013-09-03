#ifndef PARTITION_H
#define PARTITION_H

#include <complex>
typedef std::complex<long double> dcomplex;

void        neighbors(char *, int **);
int         numbp(int, int, int *);
short       bn(char);
void        initializeCanBasePair(int **);
void        translateToIntSequence(char *, short *);
void        initializeBasePairCounts(int **, int *, int );
long double hairpinloop(int, int, int, short, short, char *);
long double interiorloop(int, int, int, int, int, int, short, short, short, short);
void        solveSystem(dcomplex *, dcomplex *, char *, int **, int, int, int, int);
int         jPairedTo(int, int, int *);
int         jPairedIn(int, int, int *);
void        populateRemainingRoots(dcomplex *, dcomplex *, int, int);
void        populateMatrices(dcomplex *, dcomplex *, int, int);
void        populateZMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void        flushMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void        calculateEnergies(char *, char *, short *, int **, int **, int ***, int, double, int **, int **, int **, long double **, long double ***, long double **, long double ***, long double **, long double **, long double **);
void        evaluateZ(int, dcomplex **, dcomplex **, dcomplex **, dcomplex **, dcomplex *, dcomplex *, char *, char *, short *, int **, int **, int ***, int, int, int, int, double, dcomplex *, int **, int **, int **, long double **, long double ***, long double **, long double ***, long double **, long double **, long double **);
#endif
