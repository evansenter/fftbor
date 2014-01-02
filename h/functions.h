#ifndef FFTBOR2D_FUNCTIONS_H
#define FFTBOR2D_FUNCTIONS_H

#include "data_structures.h"
#include "mfpt/data_structures.h"

FFTBOR2D_DATA fftbor2d_from_params(FFTBOR2D_PARAMS&);
void   precalculate_energies(FFTBOR2D_DATA&);
void   evaluate_recursions_in_parallel(const FFTBOR2D_PARAMS, FFTBOR2D_DATA&, FFTBOR2D_THREADED_DATA*);
void   evaluate_recursions(int, FFTBOR2D_DATA&, FFTBOR2D_THREADED_DATA&);
void   populate_remaining_roots(FFTBOR2D_DATA&);
void   solve_system(const FFTBOR2D_PARAMS, FFTBOR2D_DATA&);
void   print_output(const FFTBOR2D_PARAMS, const FFTBOR2D_DATA);
int    j_paired_in(int, int, const int*);
void   populate_matrices(dcomplex*, int);
void   flush_matrices(dcomplex**, dcomplex**, dcomplex**, dcomplex**, int);
int*   get_bp_list(char*);
int    num_bp(int, int, const int*);
short  rna_int_code(char);
void   initialize_can_base_pair_matrix(int**);
void   translate_to_int_sequence(char*, short*);
void   initialize_base_pair_count_matrix(int**, const int*, int);
double hairpin_loop_energy(FFTBOR2D_DATA&, int, int, int, short, short, char*);
double interior_loop_energy(FFTBOR2D_DATA&, int, int, int, int, int, int, short, short, short, short);

extern "C" {
  int* get_iindx(unsigned int seq_length);
  unsigned int* maximumMatchingConstraint(const char* sequence, short* vienna_bp);
}

#endif
