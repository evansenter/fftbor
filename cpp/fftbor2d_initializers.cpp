#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits>
#include "functions.h"
#include "initializers.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define DELTA_2D(expression1, expression2, n) ((expression1) * (n) + (expression2))
// #define OPENMP_DEBUG  1
// #define SINGLE_THREAD 1

extern double temperature;

FFTBOR2D_DATA init_fftbor2d_data(FFTBOR2D_PARAMS& parameters) {
  int i, j;
  short** vienna_bp;
  int* index;
  unsigned int* max_bp_1;
  unsigned int* max_bp_2;
  int minimal_row_length;
  FFTBOR2D_DATA data = {
    NULL, // sequence
    NULL, // structure_1
    NULL, // structure_2
    0,    // seq_length
    0,    // RT
    NULL, // vienna_params
    NULL, // int_bp
    NULL, // precision_format
    NULL, // int_sequence
    NULL, // can_base_pair
    NULL, // num_base_pairs
    0,    // input_str_dist
    0,    // row_length
    0,    // run_length
    0,    // num_roots
    NULL, // solutions
    0,    // partition_function
    NULL, // roots_of_unity
    NULL, // probabilities
    NULL, // non_zero_indices
    0,    // non_zero_count
    NULL, // delta_table
    NULL, // j_paired_to_0
    NULL, // j_paired_to_1
    NULL, // EZ
    NULL, // EH
    NULL, // EHM
    NULL, // EMA
    NULL, // EMB
    NULL, // EIL
    NULL  // EM1
  };
  data.sequence    = parameters.sequence;
  data.structure_1 = parameters.structure_1;
  data.structure_2 = parameters.structure_2;
  data.seq_length  = parameters.seq_length;
  data.RT = 0.0019872370936902486 * (temperature + K0) * 100; // 0.01 * (kcal K) / mol
  read_parameter_file(parameters.energy_file);
  data.vienna_params = scale_parameters();
  data.vienna_params->model_details.special_hp = 1;
  data.int_bp    = (int**)calloc(2, sizeof(int*));
  data.int_bp[0] = get_bp_list(data.structure_1);
  data.int_bp[1] = get_bp_list(data.structure_2);
  data.precision_format = (char*)calloc(32, sizeof(char));
  sprintf(data.precision_format, "%%+.0%df", parameters.precision ? (int)floor(log(pow(2., parameters.precision)) / log(10.)) : std::numeric_limits<double>::digits);
  data.int_sequence = (short*)calloc(data.seq_length + 1, sizeof(short));
  translate_to_int_sequence(data.sequence, data.int_sequence);
  data.can_base_pair = (int**)calloc(5, sizeof(int*));

  for (i = 0; i < 5; ++i) {
    data.can_base_pair[i] = (int*)calloc(5, sizeof(int));
  }

  initialize_can_base_pair_matrix(data.can_base_pair);
  data.num_base_pairs    = (int***)calloc(2, sizeof(int**));
  data.num_base_pairs[0] = (int**)calloc(data.seq_length + 1, sizeof(int*));
  data.num_base_pairs[1] = (int**)calloc(data.seq_length + 1, sizeof(int*));

  for (i = 1; i <= data.seq_length; ++i) {
    data.num_base_pairs[0][i] = (int*)calloc(data.seq_length + 1, sizeof(int));
    data.num_base_pairs[1][i] = (int*)calloc(data.seq_length + 1, sizeof(int));
  }

  initialize_base_pair_count_matrix(data.num_base_pairs[0], data.int_bp[0], data.seq_length);
  initialize_base_pair_count_matrix(data.num_base_pairs[1], data.int_bp[1], data.seq_length);

  for (i = 1; i <= data.seq_length; ++i) {
    data.input_str_dist += (data.int_bp[0][i] > i && data.int_bp[0][i] != data.int_bp[1][i] ? 1 : 0);
    data.input_str_dist += (data.int_bp[1][i] > i && data.int_bp[1][i] != data.int_bp[0][i] ? 1 : 0);
  }

  // Secondary structure data structure in the slightly different format that Vienna uses (1-indexed short array with 0 as unpaired sentinel value).
  vienna_bp = (short**)calloc(2, sizeof(short*));

  for (i = 0; i < 2; ++i) {
    vienna_bp[i]    = (short*)calloc(data.seq_length + 1, sizeof(short));;
    vienna_bp[i][0] = data.seq_length;

    for (j = 1; j <= data.seq_length; ++j) {
      vienna_bp[i][j] = data.int_bp[i][j] > 0 ? data.int_bp[i][j] : 0;
    }
  }

  // Index for moving in quadratic distancy dimensions (from Vienna 2.1.2)
  index = get_iindx((unsigned)data.seq_length);
  // Maximally saturated structure constrained with the input structures.
  max_bp_1 = maximumMatchingConstraint(data.sequence, vienna_bp[0]);
  max_bp_2 = maximumMatchingConstraint(data.sequence, vienna_bp[1]);
  // Minimize the row size to the max BP distance.
  minimal_row_length = MAX2(data.int_bp[0][0] + max_bp_1[index[1] - data.seq_length], data.int_bp[1][0] + max_bp_2[index[1] - data.seq_length]);
  // Note: row_length = (least even number >= minimal_row_length) + 1 (for a seq. with minimal_row_length = 9, row_length = (0..10).length = 11)
  data.row_length = (minimal_row_length % 2 ? minimal_row_length + 1 : minimal_row_length) + 1;
  // Note: run_length = (least number div. 4 >= row_length ^ 2) / 2 (for a seq. with row_length = 11, run_length = (11 ^ 2 + 3) / 2 = 62)
  data.run_length = ((int)pow((double)data.row_length, 2) + ((int)pow((double)data.row_length, 2) % 4)) / 2;
  data.num_roots  = data.run_length * 2;
  // Initialize convenience tables for storing solutions, roots of unity, probabilities and non-zero indices.
  data.solutions        = (dcomplex*)calloc(data.run_length + 1, sizeof(dcomplex));
  data.roots_of_unity   = (dcomplex*)calloc(data.num_roots, sizeof(dcomplex));
  data.probabilities    = (double*)calloc(2 * data.run_length + 1, sizeof(double));
  data.non_zero_indices = (int*)calloc((int)pow(data.row_length, 2.) + 1, sizeof(int));
  // Populate convenience tables.
  populate_matrices(data.roots_of_unity, data.num_roots);
  // Create convenience table for looking up 1D indexing of (k, l) coordinates.
  data.delta_table = (int**)calloc(data.seq_length + 1, sizeof(int*));

  for (i = 0; i <= data.seq_length; ++i) {
    data.delta_table[i] = (int*)calloc(data.seq_length + 1, sizeof(int));

    for (j = 0; j <= data.seq_length; ++j) {
      data.delta_table[i][j] = DELTA_2D(i, j, data.row_length);
    }
  }

  // Create convenience table for boolean (i, j paired?) values.
  data.j_paired_to_0 = (int**)calloc(data.seq_length + 1, sizeof(int*));
  data.j_paired_to_1 = (int**)calloc(data.seq_length + 1, sizeof(int*));

  for (i = 0; i <= data.seq_length; ++i) {
    data.j_paired_to_0[i] = (int*)calloc(data.seq_length + 1, sizeof(int));
    data.j_paired_to_1[i] = (int*)calloc(data.seq_length + 1, sizeof(int));

    for (j = 0; j <= data.seq_length; ++j) {
      data.j_paired_to_0[i][j] = j_paired_to(i, j, data.int_bp[0]);
      data.j_paired_to_1[i][j] = j_paired_to(i, j, data.int_bp[1]);
    }
  }

  // Initialize tables for precalculating energies.
  data.EZ  = (double**)calloc(data.seq_length + 1, sizeof(double*));
  data.EH  = (double**)calloc(data.seq_length + 1, sizeof(double*));
  data.EHM = (double**)calloc(data.seq_length + 1, sizeof(double*));
  data.EMA = (double**)calloc(data.seq_length + 1, sizeof(double*));
  data.EMB = (double**)calloc(data.seq_length + 1, sizeof(double*));
  data.EIL = (double***)calloc(data.seq_length + 1, sizeof(double**));;
  data.EM1 = (double***)calloc(data.seq_length + 1, sizeof(double**));;

  for (i = 0; i <= data.seq_length; ++i) {
    data.EZ[i]  = (double*)calloc(data.seq_length + 1, sizeof(double));
    data.EH[i]  = (double*)calloc(data.seq_length + 1, sizeof(double));
    data.EHM[i] = (double*)calloc(data.seq_length + 1, sizeof(double));
    data.EMA[i] = (double*)calloc(data.seq_length + 1, sizeof(double));
    data.EMB[i] = (double*)calloc(data.seq_length + 1, sizeof(double));
    data.EIL[i] = (double**)calloc(data.seq_length + 1, sizeof(double*));
    data.EM1[i] = (double**)calloc(data.seq_length + 1, sizeof(double*));

    for (j = 0; j <= data.seq_length; ++j) {
      data.EM1[i][j] = (double*)calloc(data.seq_length + 1, sizeof(double));
      // Multiplied by a "twiddle" factor.
      data.EIL[i][j] = (double*)calloc(6 * (data.seq_length + 1), sizeof(double));
    }
  }

  return data;
}

void free_fftbor2d_data(FFTBOR2D_DATA& data) {
  free(data.sequence);
  free(data.structure_1);
  free(data.structure_2);
  free(data.vienna_params);
  free(data.int_bp);
  free(data.precision_format);
  free(data.int_sequence);
  free(data.can_base_pair);
  free(data.num_base_pairs);
  free(data.solutions);
  free(data.roots_of_unity);
  free(data.probabilities);
  free(data.non_zero_indices);
  free(data.delta_table);
  free(data.j_paired_to_0);
  free(data.j_paired_to_1);
  free(data.EZ);
  free(data.EH);
  free(data.EHM);
  free(data.EMA);
  free(data.EMB);
  free(data.EIL);
  free(data.EM1);
}

void print_fftbor2d_data(FFTBOR2D_DATA& data) {
  printf("FFTBOR2D_DATA:\n");
  printf("sequence\t\t%s\n",         data.sequence    == NULL ? "*missing*" : data.sequence);
  printf("structure_1\t\t%s\n",      data.structure_1 == NULL ? "*missing*" : data.structure_1);
  printf("structure_2\t\t%s\n",      data.structure_2 == NULL ? "*missing*" : data.structure_2);
  printf("seq_length\t\t%d\n",       data.seq_length);
  printf("RT\t\t\t%f\n",             data.RT);
  printf("precision_format\t%s\n",   data.precision_format == NULL ? "*missing*" : data.precision_format);
  printf("input_str_dist\t\t%d\n",   data.input_str_dist);
  printf("row_length\t\t%d\n",       data.row_length);
  printf("run_length\t\t%d\n",       data.run_length);
  printf("num_roots\t\t%d\n",        data.num_roots);
  printf("partition_function\t%f\n", data.partition_function);
  printf("non_zero_count\t\t%d\n",   data.non_zero_count);
  printf("\n");
}

FFTBOR2D_THREADED_DATA* init_fftbor2d_threaded_data(FFTBOR2D_PARAMS& parameters, FFTBOR2D_DATA& data) {
  int i, j;
  FFTBOR2D_THREADED_DATA* threaded_data;
  #if defined(_OPENMP) && defined(OPENMP_DEBUG)
  printf("Max threads possible: %d\n", omp_get_max_threads());
  #ifdef SINGLE_THREAD
  parameters.max_threads = 1;
  #endif
  printf("Setting number of threads: %d\n", parameters.max_threads);
  #endif
  #ifdef _OPENMP
  // Set number of threads for OpenMP
  omp_set_num_threads(parameters.max_threads);
  #endif
  threaded_data = (FFTBOR2D_THREADED_DATA*)calloc(parameters.max_threads, sizeof(FFTBOR2D_THREADED_DATA));

  for (i = 0; i < parameters.max_threads; ++i) {
    threaded_data[i] = {
      NULL, // Z
      NULL, // ZB
      NULL, // ZM
      NULL, // ZM1
      NULL  // root_to_power
    };
    // Initialize matricies for dynamic programming.
    threaded_data[i].Z   = (dcomplex**)calloc(data.seq_length + 1, sizeof(dcomplex*));
    threaded_data[i].ZB  = (dcomplex**)calloc(data.seq_length + 1, sizeof(dcomplex*));
    threaded_data[i].ZM  = (dcomplex**)calloc(data.seq_length + 1, sizeof(dcomplex*));
    threaded_data[i].ZM1 = (dcomplex**)calloc(data.seq_length + 1, sizeof(dcomplex*));

    for (j = 0; j <= data.seq_length; ++j) {
      threaded_data[i].Z[j]   = (dcomplex*)calloc(data.seq_length + 1, sizeof(dcomplex));
      threaded_data[i].ZB[j]  = (dcomplex*)calloc(data.seq_length + 1, sizeof(dcomplex));
      threaded_data[i].ZM[j]  = (dcomplex*)calloc(data.seq_length + 1, sizeof(dcomplex));
      threaded_data[i].ZM1[j] = (dcomplex*)calloc(data.seq_length + 1, sizeof(dcomplex));
    }

    threaded_data[i].root_to_power = (dcomplex*)calloc(data.num_roots, sizeof(dcomplex));
  }

  return threaded_data;
}

void free_fftbor2d_threaded_data(FFTBOR2D_THREADED_DATA* threaded_data, int max_threads) {
  int i;

  for (i = 0; i < max_threads; ++i) {
    free(threaded_data[i].Z);
    free(threaded_data[i].ZB);
    free(threaded_data[i].ZM);
    free(threaded_data[i].ZM1);
    free(threaded_data[i].root_to_power);
  }

  free(threaded_data);
}

inline int j_paired_to(int i, int j, int* base_pairs) {
  return base_pairs[i] == j ? -1 : 1;
}
