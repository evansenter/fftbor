#ifndef PARTITION_H
#define PARTITION_H

#include <complex>
typedef std::complex<double> dcomplex;

void   neighbors(char *, int **);
int    numbp(int, int, int *);
short  bn(char);
void   initializeCanBasePair(int **);
void   translateToIntSequence(char *, short *);
void   initializeBasePairCounts(int **, int *, int );
double hairpinloop(int, int, int, short, short, char *);
double interiorloop(int, int, int, int, int, int, short, short, short, short);
void   solveSystem(double *, dcomplex *, int, int, int, int *, int&, double&, int&, char *);
void   printOutput(double *, int, int, int, int *, int&, double&, int&, char *);
void   computeTransitionMatrix(int *, int, double *, double **);
double computeMFPT(int *, int, double **, int, int);
int    jPairedTo(int, int, int *);
int    jPairedIn(int, int, int *);
void   populateRemainingRoots(dcomplex *, dcomplex *, int, int);
void   populateMatrices(dcomplex *, int);
void   populateZMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void   flushMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void   calculateEnergies(char *, short *, int **, int, double, double **, double ***, double **, double ***, double **, double **, double **);
void   evaluateZ(int, dcomplex **, dcomplex **, dcomplex **, dcomplex **, dcomplex *, dcomplex *, char *, short *, int **, int **, int ***, int, int, int, dcomplex *, int **, int **, int **,double **, double ***, double **, double ***, double **, double **, double **);
void   inverse(double*, int);
char*  findEnergyFile();
#endif
