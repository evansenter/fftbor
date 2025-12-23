#ifndef DELTA_H
#define DELTA_H

#include <complex>
#include <vector>
#include "memory_types.h"

// Use the dcomplex type from fftbor namespace (defined in memory_types.h)
// Provide a global alias for backward compatibility
using dcomplex = fftbor::dcomplex;

// Main entry point
void neighbours(const char* input_sequence, const int* bp_list);

// Base pair utilities
int num_bp(int i, int j, const int* bp_list);
int encode_base(char base);

// Initialization functions
void initialize_can_base_pair(fftbor::IntMatrix2D& can_base_pair);
void translate_to_int_sequence(const char* sequence, int* int_sequence);
void initialize_base_pair_counts(fftbor::IntMatrix2D& num_base_pairs, const int* bp_list, int n);

// Energy calculations
double hairpin_loop(int i, int j, int type, short si1, short sj1, const char* string);
double interior_loop(int i, int j, int k, int l, int type, int type_2, short si1, short sq1, short sj1, short sp1);

// Partition function solver
void solve_system(fftbor::ComplexMatrix3D& solutions, const char* sequence, const int* structure, int sequence_length, int run_length);

// Base pair query functions
int j_paired_to(int i, int j, const int* base_pairs);
int j_paired_in(int i, int j, const int* base_pairs);

// Matrix operations
void populate_remaining_roots(fftbor::ComplexMatrix3D& solutions, int sequence_length, int run_length, int last_root);

void populate_matrices(fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                       fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                       fftbor::ComplexMatrix3D& solutions, std::vector<dcomplex>& roots_of_unity,
                       int sequence_length, int run_length);

void flush_matrices(fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                    fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                    int sequence_length);

void evaluate_z(int root,
                fftbor::ComplexMatrix2D& Z, fftbor::ComplexMatrix2D& ZB,
                fftbor::ComplexMatrix2D& ZM, fftbor::ComplexMatrix2D& ZM1,
                fftbor::ComplexMatrix3D& solutions, const std::vector<dcomplex>& roots_of_unity,
                const char* input_sequence, const char* sequence,
                const int* int_sequence, const int* bp_list,
                const fftbor::IntMatrix2D& can_base_pair, const fftbor::IntMatrix2D& num_base_pairs,
                int sequence_length, int run_length, double RT);

#endif
