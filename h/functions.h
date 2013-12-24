#ifndef FFTBOR2D_FUNCTIONS_H
#define FFTBOR2D_FUNCTIONS_H

#include "mfpt/data_structures.h"
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
void   solveSystem(double *, dcomplex *, int, int, int, int *, int&, double&, char *);
void   printOutput(double *, int, int, int, int *, int&, double&, char *);
void   calculateKinetics(int *, int, double *, int, char *);
void   convert_fftbor2d_energy_grid_to_klp_matrix(int *, double *, int, KLP_MATRIX);
void   populationProportion(int *, int, double *, int, char *);
double* convert_fftbor2d_energy_grid_to_transition_rate_matrix(int *, int, double *);
int    jPairedTo(int, int, int *);
int    jPairedIn(int, int, int *);
void   populateRemainingRoots(dcomplex *, dcomplex *, int, int);
void   populateMatrices(dcomplex *, int);
void   populateZMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void   flushMatrices(dcomplex **, dcomplex **, dcomplex **, dcomplex **, int);
void   calculateEnergies(char *, short *, int **, int, double, double **, double ***, double **, double ***, double **, double **, double **);
void   evaluateZ(int, dcomplex **, dcomplex **, dcomplex **, dcomplex **, dcomplex *, dcomplex *, char *, short *, int **, int **, int ***, int, int, int, dcomplex *, int **, int **, int **,double **, double ***, double **, double ***, double **, double **, double **);
char*  findEnergyFile();
#endif
