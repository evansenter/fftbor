#include "params.h"

#ifndef DELTA_H
#define DELTA_H

/* To avoid computing MFE structures if this is not wanted. */
//#define COMPUTEMFE


int DELTA;
int PF;
int DANGLE;
int NUMBER;
char *ENERGY;
int STRUCTURE;
int PARTITION;
int STOP;
int N;
extern double temperature;

/* Returns 1 if the two input characters are complementary bases */
int basepair(char , char );

void sub_structure(int , int , int *, int *);
void neighbours(char *,int *);
void pf(char *);
int basepaired_to(int ,int *);
int bp_diff(int *, int, int , int);
void print_bps(int *);

paramT *P;

#endif
