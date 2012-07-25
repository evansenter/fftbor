#ifndef MISC_H
#define MISC_H

#include <stdbool.h>

#define MAX_SEQ_LENGTH 610

void *xcalloc(size_t, size_t);
int *getBasePairList(char *);
int min2(int, int);
int max2(int, int);

#endif
