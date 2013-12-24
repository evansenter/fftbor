#ifndef MISC_H
#define MISC_H

#include <time.h>
#include <sys/time.h>
#include <stdbool.h>

#define MAX_SEQ_LENGTH 610
#define TRIEQUALS(x, y, z) ((x == y) && (y == z)) /* Transitivity (BOOM) */

void *xcalloc(size_t, size_t);
int *getBasePairList(char *);
int min2(int, int);
int max2(int, int);

#endif
