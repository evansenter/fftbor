#ifndef MISC_H
#define MISC_H

#define DEBUG 0

/* Constant definitions */
#define THRESHOLD  3
#define INF 1000000
#define AUenergy  -2
#define GUenergy  -1
#define GCenergy  -3
#define NUM        200000  /* Upper bound on number of Monte Carlo steps */
//#define NUMRUNS    100  /* Upper bound on number of runs of entire folding
//			 * process from empty sec structure to native state */
#define NN 610 /* Maximum length of RNA */

#include <stdbool.h>

void *xcalloc(size_t , size_t );
int *getBasePairList(char *);
int mindouble(double, double *);
int min2(int , int);
int max2(int , int);

#endif
