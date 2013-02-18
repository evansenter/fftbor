#ifndef MISC_H
#define MISC_H

#include <time.h>
#include <sys/time.h>
#include <stdbool.h>

#define TIMING(start, stop, task) printf("Time in ms for %s: %.2f\n", task, (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0));
#define MAX_SEQ_LENGTH 610
#define TRIEQUALS(x, y, z) ((x == y) && (y == z)) /* Transitivity (BOOM) */

void *xcalloc(size_t, size_t);
int *getBasePairList(char *);
int min2(int, int);
int max2(int, int);

#endif
