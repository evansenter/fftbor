#ifndef MISC_H
#define MISC_H

#include <cstddef>  // for size_t
#include "memory_types.h"

constexpr int MAX_SEQ_LENGTH = 610;

void *xcalloc(size_t n, size_t s);
fftbor::BasePairListPtr getBasePairList(const char* secStr);
int min2(int m, int n);
int max2(int m, int n);

#endif
