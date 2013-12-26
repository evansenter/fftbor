#ifndef FFTBOR2D_DATA_STRUCTURES_H
#define FFTBOR2D_DATA_STRUCTURES_H

#include "vienna/data_structures.h"
#include <complex>
typedef std::complex<double> dcomplex;

typedef struct {
  char* sequence;
  char* structure_1;
  char* structure_2;
  char* energy_file;
  int   seq_length;
  int   precision;
  int   max_threads;
  char  format;
  short verbose;
} FFTBOR2D_PARAMS;

typedef struct { // Variables are sorted by the order they get instantiated.
  char*     sequence;
  char*     structure_1;
  char*     structure_2;
  int       seq_length;
  double    RT;
  paramT*   vienna_params;
  int**     int_bp;
  char*     precision_format;
  short*    int_sequence;
  int**     can_base_pair;
  int***    num_base_pairs;
  int       input_str_dist;
  int       row_length;
  int       run_length;
  int       num_roots;
  dcomplex* solutions;
  double    partition_function;
  dcomplex* roots_of_unity;
  double*   probabilities;
  int*      non_zero_indices;
  int       non_zero_count;
  int**     delta_table;
  int**     j_paired_to_0;
  int**     j_paired_to_1;
  double**  EZ;
  double**  EH;
  double**  EHM;
  double**  EMA;
  double**  EMB;
  double*** EIL;
  double*** EM1;
} FFTBOR2D_DATA;

typedef struct { // Variables are sorted by the order they get instantiated.
  dcomplex** Z;
  dcomplex** ZB;
  dcomplex** ZM;
  dcomplex** ZM1;
  dcomplex*  root_to_power;
} FFTBOR2D_THREADED_DATA;

#endif
