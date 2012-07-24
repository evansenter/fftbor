#ifndef DELTA_H
#define DELTA_H

#include "params.h"
#include <complex>

typedef std::complex<double> dcomplex;

void   neighbours(char *, int *);
void   pf(char *);
int    numbp(int, int, int *);
int    bn(char);
void   initialize_PN(int[][5]);
void   translateseq(char *, int *);
void   initialize_NumBP(int **, int *, int );
double hairpinloop(int, int, int, int *, char *);
double interiorloop(int, int, int, int, int, int, int *);
double multiloop_closing(int, int, int, int, int, int, int *);
void   solveSystem(dcomplex **, double *, double, int);
int    jPairedTo(int, int, int *);
int    jPairedIn(int, int, int *);

#endif