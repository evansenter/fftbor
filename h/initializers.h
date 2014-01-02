#ifndef FFTBOR2D_INITIALIZERS_H
#define FFTBOR2D_INITIALIZERS_H

#define MIN2(A, B) ((A) < (B) ? (A) : (B))
#define MAX2(A, B) ((A) > (B) ? (A) : (B))

#include "data_structures.h"

FFTBOR2D_DATA init_fftbor2d_data(const FFTBOR2D_PARAMS);
void free_fftbor2d_data(FFTBOR2D_DATA&);
int minimum_row_length(const FFTBOR2D_DATA);
void print_fftbor2d_data(const FFTBOR2D_DATA);
FFTBOR2D_THREADED_DATA* init_fftbor2d_threaded_data(FFTBOR2D_PARAMS&, const FFTBOR2D_DATA);
void free_fftbor2d_threaded_data(FFTBOR2D_THREADED_DATA*, int);
int j_paired_to(int, int, int*);

extern "C" {
  void read_parameter_file(const char*);
  paramT* scale_parameters(void);
}

#endif
