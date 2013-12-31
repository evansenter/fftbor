#ifndef FFTBOR2D_PARAMS_H
#define FFTBOR2D_PARAMS_H

#include "data_structures.h"

FFTBOR2D_PARAMS init_fftbor2d_params();
void free_fftbor2d_params(FFTBOR2D_PARAMS);
FFTBOR2D_PARAMS parse_fftbor2d_args(int, char**);
void parse_fftbor2d_sequence_data(int, char**, int, FFTBOR2D_PARAMS&);
char* find_energy_file(char*);
int fftbor2d_error_handling(FFTBOR2D_PARAMS&);
void debug_fftbor2d_parameters(FFTBOR2D_PARAMS&);
void fftbor2d_usage();

#endif
