#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "params.h"
#include "initializers.h"
#include "functions.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MIN_PAIR_DIST TURN
#define MAX_INTERIOR_DIST MAXLOOP
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define FFTW_REAL 0
#define FFTW_IMAG 1
#define COMPLEX_CONJ(complexNumber) (dcomplex((complexNumber).real(), -(complexNumber).imag()))
#define ROOT_POW(i, pow, n) (data.roots_of_unity[((i) * (pow)) % (n)])
#define PRINT_COMPLEX(i, complex) printf("%d: %+f %+fi\n", i, complex[i].real(), complex[i].imag())

// #define SILENCE_OUTPUT  1
// #define FFTBOR_DEBUG    1
// #define TWIDDLE_DEBUG   1
// #define MEASURE_TWIDDLE 1
// #define STRUCTURE_COUNT 1
#define DO_WORK

void precalculate_energies(FFTBOR2D_DATA& data) {
  int i, j, k, l, d, position;
  #if MEASURE_TWIDDLE
  double max_twiddle = 0;
  #endif
  
  for (i = 1; i <= data.seq_length; ++i) {
    for (d = MIN_PAIR_DIST + 1; d <= data.seq_length - i; ++d) {
      j = i + d;
      
      if (data.can_base_pair[data.int_sequence[i]][data.int_sequence[j]]) {
        // ****************************************************************************
        // Solve ZB
        // ****************************************************************************
        // In a hairpin, [i + 1, j - 1] unpaired.
        data.EH[i][j] = exp(-hairpin_loop_energy(data, i, j, data.can_base_pair[data.int_sequence[i]][data.int_sequence[j]], data.int_sequence[i + 1], data.int_sequence[j - 1], data.sequence) / data.RT);
        position = 0;
        
        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < MIN2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          // There can't be more than 30 unpaired bases between i and k, and there must be room between k and j for l
          for (l = MAX2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) {
            // l needs to at least have room to pair with k, and there can be at most 30 unpaired bases between (i, k) + (l, j), with l < j
            if (data.can_base_pair[data.int_sequence[k]][data.int_sequence[l]]) {
              #ifdef TWIDDLE_DEBUG
            
              if (position >= 6 * (data.seq_length + 1)) {
                fprintf(stderr, "Trying to access non-existant memory in EIL at index %d (%f).\n", position, position / (data.seq_length + 1.));
              }
              
              #endif
              #ifdef MEASURE_TWIDDLE
              max_twiddle = position / (data.seq_length + 1.) > max_twiddle ? position / (data.seq_length + 1.) : max_twiddle;
              #endif
              // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
              data.EIL[i][j][position] = exp(-interior_loop_energy(data, i, j, k, l, data.can_base_pair[data.int_sequence[i]][data.int_sequence[j]], data.can_base_pair[data.int_sequence[l]][data.int_sequence[k]], data.int_sequence[i + 1], data.int_sequence[l + 1], data.int_sequence[j - 1], data.int_sequence[k - 1]) / data.RT);
              position++;
            }
          }
        }
        
        // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1)
        data.EHM[i][j] = exp(-(data.vienna_params->MLclosing + data.vienna_params->MLintern[data.can_base_pair[data.int_sequence[i]][data.int_sequence[j]]]) / data.RT);
      }
      
      // ****************************************************************************
      // Solve ZM1
      // ****************************************************************************
      for (k = i + MIN_PAIR_DIST + 1; k <= j; ++k) {
        // k is the closing base pairing with i of a single component within the range [i, j]
        if (data.can_base_pair[data.int_sequence[i]][data.int_sequence[k]]) {
          data.EM1[i][j][k] = exp(-(data.vienna_params->MLbase * (j - k) + data.vienna_params->MLintern[data.can_base_pair[data.int_sequence[i]][data.int_sequence[k]]]) / data.RT);
        }
      }
      
      // ****************************************************************************
      // Solve ZM
      // ****************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        // Only one stem.
        data.EMA[i][k] = exp(-(data.vienna_params->MLbase * (k - i)) / data.RT);
        
        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
          data.EMB[j][k] = exp(-(data.vienna_params->MLintern[data.can_base_pair[data.int_sequence[k]][data.int_sequence[j]]]) / data.RT);
        }
      }
      
      // **************************************************************************
      // Solve Z
      // **************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        // (k, j) is the rightmost base pair in (i, j)
        if (data.can_base_pair[data.int_sequence[k]][data.int_sequence[j]]) {
          data.EZ[j][k] = exp(-(data.can_base_pair[data.int_sequence[k]][data.int_sequence[j]] > 2 ? data.vienna_params->TerminalAU : 0) / data.RT);
        }
      }
    }
  }
  
  #ifdef MEASURE_TWIDDLE
  printf("Max twiddle distance seen: %f\n", max_twiddle);
  #endif
}

void evaluate_recursions_in_parallel(FFTBOR2D_PARAMS& parameters, FFTBOR2D_DATA& data, FFTBOR2D_THREADED_DATA* threaded_data) {
  int i, thread_id;
  
  // Start main recursions (i <= data.run_length / 2 is an optimization leveraging complex conjugates).
  #pragma omp parallel for private(i, thread_id) shared(data, threaded_data) default(none) num_threads(parameters.max_threads)

  for (i = 0; i <= data.run_length / 2; ++i) {
    #ifdef _OPENMP
    thread_id = omp_get_thread_num();
    #else
    thread_id = 0;
    #endif
    
    evaluate_recursions(i, data, threaded_data[thread_id]);
  }
}

void evaluate_recursions(int root, FFTBOR2D_DATA& data, FFTBOR2D_THREADED_DATA& threaded_data) {
  int i, j, k, l, d, delta, position;
  double energy;
  flush_matrices(threaded_data.Z, threaded_data.ZB, threaded_data.ZM, threaded_data.ZM1, data.seq_length);
  
  for (i = 0; i < data.num_roots; ++i) {
    threaded_data.root_to_power[i] = ROOT_POW(root, i, data.num_roots);
  }
  
  for (d = MIN_PAIR_DIST + 1; d < data.seq_length; ++d) {
    for (i = 1; i <= data.seq_length - d; ++i) {
      j = i + d;
      
      if (data.can_base_pair[data.int_sequence[i]][data.int_sequence[j]]) {
        // ****************************************************************************
        // Solve ZB
        // ****************************************************************************
        #ifdef DO_WORK
        // In a hairpin, [i + 1, j - 1] unpaired.
        delta  = data.delta_table[data.num_base_pairs[0][i][j] + data.j_paired_to_0[i][j]][data.num_base_pairs[1][i][j] + data.j_paired_to_1[i][j]];
        threaded_data.ZB[i][j] += threaded_data.root_to_power[delta] * data.EH[i][j];
        #ifdef STRUCTURE_COUNT
        threaded_data.ZB[j][i] += 1;
        #endif
        #endif
        position = 0;
        
        // Interior loop / bulge / stack / multiloop.
        for (k = i + 1; k < MIN2(i + 30, j - MIN_PAIR_DIST - 2) + 1; ++k) {
          // There can't be more than 30 unpaired bases between i and k, and there must be room between k and j for l
          for (l = MAX2(k + MIN_PAIR_DIST + 1, j - (MAX_INTERIOR_DIST - (k - i))); l < j; ++l) {
            // l needs to at least have room to pair with k, and there can be at most 30 unpaired bases between (i, k) + (l, j), with l < j
            if (data.can_base_pair[data.int_sequence[k]][data.int_sequence[l]]) {
              #ifdef DO_WORK
              // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
              delta  = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][k][l] + data.j_paired_to_0[i][j]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][k][l] + data.j_paired_to_1[i][j]];
              threaded_data.ZB[i][j] += threaded_data.ZB[k][l] * threaded_data.root_to_power[delta] * data.EIL[i][j][position];
              position++;
              #ifdef STRUCTURE_COUNT
              threaded_data.ZB[j][i] += threaded_data.ZB[l][k];
              #endif
              #endif
            }
          }
        }
        
        energy = data.EHM[i][j];
        
        for (k = i + MIN_PAIR_DIST + 3; k < j - MIN_PAIR_DIST - 1; ++k) {
          #ifdef DO_WORK
          // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1)
          delta  = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][i + 1][k - 1] - data.num_base_pairs[0][k][j - 1] + data.j_paired_to_0[i][j]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][i + 1][k - 1] - data.num_base_pairs[1][k][j - 1] + data.j_paired_to_1[i][j]];
          threaded_data.ZB[i][j] += threaded_data.ZM[i + 1][k - 1] * threaded_data.ZM1[k][j - 1] * threaded_data.root_to_power[delta] * energy;
          #ifdef STRUCTURE_COUNT
          threaded_data.ZB[j][i] += threaded_data.ZM[k - 1][i + 1] * threaded_data.ZM1[j - 1][k];
          #endif
          #endif
        }
      }
      
      // ****************************************************************************
      // Solve ZM1
      // ****************************************************************************
      for (k = i + MIN_PAIR_DIST + 1; k <= j; ++k) {
        // k is the closing base pairing with i of a single component within the range [i, j]
        if (data.can_base_pair[data.int_sequence[i]][data.int_sequence[k]]) {
          #ifdef DO_WORK
          delta  = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][i][k]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][i][k]];
          threaded_data.ZM1[i][j] += threaded_data.ZB[i][k] * data.EM1[i][j][k];
          #ifdef STRUCTURE_COUNT
          threaded_data.ZM1[j][i] += threaded_data.ZB[k][i];
          #endif
          #endif
        }
      }
      
      // ****************************************************************************
      // Solve ZM
      // ****************************************************************************
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        #ifdef DO_WORK
        // Only one stem.
        delta  = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][k][j]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][k][j]];
        threaded_data.ZM[i][j] += threaded_data.ZM1[k][j] * threaded_data.root_to_power[delta] * data.EMA[i][k];
        #ifdef STRUCTURE_COUNT
        threaded_data.ZM[j][i] += threaded_data.ZM1[j][k];
        #endif
        #endif
        
        // More than one stem.
        if (k > i + MIN_PAIR_DIST + 1) { // (k > i + MIN_PAIR_DIST + 1) because i can pair with k - 1
          #ifdef DO_WORK
          delta  = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][i][k - 1] - data.num_base_pairs[0][k][j]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][i][k - 1] - data.num_base_pairs[1][k][j]];
          threaded_data.ZM[i][j] += threaded_data.ZM[i][k - 1] * threaded_data.ZM1[k][j] * threaded_data.root_to_power[delta] * data.EMB[j][k];
          #ifdef STRUCTURE_COUNT
          threaded_data.ZM[j][i] += threaded_data.ZM[k - 1][i] * threaded_data.ZM1[j][k];
          #endif
          #endif
        }
      }
      
      // **************************************************************************
      // Solve Z
      // **************************************************************************
      #ifdef DO_WORK
      delta = data.delta_table[j_paired_in(i, j, data.int_bp[0])][j_paired_in(i, j, data.int_bp[1])];
      threaded_data.Z[i][j] += threaded_data.Z[i][j - 1] * threaded_data.root_to_power[delta];
      #ifdef STRUCTURE_COUNT
      threaded_data.Z[j][i] += threaded_data.Z[j - 1][i];
      #endif
      #endif
      
      for (k = i; k < j - MIN_PAIR_DIST; ++k) {
        // (k, j) is the rightmost base pair in (i, j)
        if (data.can_base_pair[data.int_sequence[k]][data.int_sequence[j]]) {
          #ifdef DO_WORK
          energy = data.EZ[j][k];
          #endif
          
          if (k == i) {
            #ifdef DO_WORK
            delta = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][k][j]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][k][j]];
            threaded_data.Z[i][j] += threaded_data.ZB[k][j] * threaded_data.root_to_power[delta] * energy;
            #ifdef STRUCTURE_COUNT
            threaded_data.Z[j][i] += threaded_data.ZB[j][k];
            #endif
            #endif
          } else {
            #ifdef DO_WORK
            delta = data.delta_table[data.num_base_pairs[0][i][j] - data.num_base_pairs[0][i][k - 1] - data.num_base_pairs[0][k][j]][data.num_base_pairs[1][i][j] - data.num_base_pairs[1][i][k - 1] - data.num_base_pairs[1][k][j]];
            threaded_data.Z[i][j] += threaded_data.Z[i][k - 1] * threaded_data.ZB[k][j] * threaded_data.root_to_power[delta] * energy;
            #ifdef STRUCTURE_COUNT
            threaded_data.Z[j][i] += threaded_data.Z[k - 1][i] * threaded_data.ZB[j][k];
            #endif
            #endif
          }
        }
      }
    }
  }
  
  data.solutions[root] = threaded_data.Z[1][data.seq_length];
  #ifdef FFTBOR_DEBUG
  printf(".");
  #endif
}

void populate_remaining_roots(FFTBOR2D_DATA& data) {
  // Optimization leveraging complex conjugate of data.solutions to the polynomial.
  int i;
  #ifdef FFTBOR_DEBUG
  printf("Solutions (before populating remaining roots):\n");
  
  for (i = 0; i < data.run_length; ++i) {
    PRINT_COMPLEX(i, data.solutions);
  }
  
  #endif
  
  for (i = data.run_length / 2 + 1; i < data.run_length; ++i) {
    #ifdef FFTBOR_DEBUG
    printf("l: %d, r %d\n", i, data.run_length - i);
    #endif
    data.solutions[i] = COMPLEX_CONJ(data.solutions[data.run_length - i]);
  }
  
  #ifdef FFTBOR_DEBUG
  printf("Solutions (after populating remaining roots):\n");
  
  for (i = 0; i < data.run_length; ++i) {
    PRINT_COMPLEX(i, data.solutions);
  }
  
  #endif
  
  if (data.input_str_dist % 2) {
    for (i = 0; i < data.run_length; ++i) {
      data.solutions[i] = COMPLEX_CONJ(data.roots_of_unity[i]) * data.solutions[i];
    }
    
    for (i = data.run_length / 2 + 1; i <= data.run_length; ++i) {
      data.solutions[i] = dcomplex(-1, 0) * data.solutions[i];
    }
    
    #ifdef FFTBOR_DEBUG
    printf("Solutions (after multiplying data.solutions by nu^{-k} (and the scalar -1 for data.solutions y_{k: k > M_{0} / 2})):\n");
    
    for (i = 0; i < data.run_length; ++i) {
      PRINT_COMPLEX(i, data.solutions);
    }
    
    #endif
  }
}

void solve_system(FFTBOR2D_PARAMS& parameters, FFTBOR2D_DATA& data) {
  int i, x, y;
  int offset              = data.input_str_dist % 2 ? 1 : 0;
  double sum              = 0;
  double normal_sum       = 0;
  fftw_complex* signal    = (fftw_complex*)calloc(data.run_length, sizeof(fftw_complex));
  fftw_complex* result    = (fftw_complex*)calloc(data.run_length, sizeof(fftw_complex));
  fftw_plan plan          = fftw_plan_dft_1d(data.run_length, signal, result, FFTW_BACKWARD, FFTW_ESTIMATE);
  data.partition_function = data.solutions[0].real();
  
  // For some reason it's much more numerically stable to set the signal via real / imag components separately.
  for (i = 0; i < data.run_length; ++i) {
    // Convert point-value data.solutions of Z(root) to 10^parameters.precision * Z(root) / Z
    signal[i][FFTW_REAL] = pow(2., (double)parameters.precision) * (data.solutions[i].real() / data.partition_function);
    signal[i][FFTW_IMAG] = pow(2., (double)parameters.precision) * (data.solutions[i].imag() / data.partition_function);
  }
  
  #ifdef FFTBOR_DEBUG
  printf("Scaling factor: %f:\n", data.solutions[0].real());
  printf("Scaled data.solutions vector for the inverse DFT:\n");
  
  for (i = 0; i < data.run_length; ++i) {
    printf("%d: %+f %+fi\n", i, signal[i][FFTW_REAL], signal[i][FFTW_IMAG]);
  }
  
  #endif
  // Calculate transform, coefficients are in fftw_complex result array.
  fftw_execute(plan);
  
  for (i = 0; i < data.run_length; ++i) {
    // Truncate to user-specified precision; if set to 0, no truncation occurs (dangerous).
    if (!parameters.precision) {
      data.solutions[i] = dcomplex(result[i][FFTW_REAL] / data.run_length, 0);
    } else {
      data.solutions[i] = dcomplex(pow(2., -parameters.precision) * static_cast<int>(result[i][FFTW_REAL] / data.run_length), 0);
    }
    
    x = (2 * i + offset) / data.row_length;
    y = (2 * i + offset) % data.row_length;
    
    // Probabilities must be > 0 and satisfy the triangle inequality.
    if (
      data.solutions[i].real() > 0 &&
      x + y >= data.input_str_dist &&
      x + data.input_str_dist >= y &&
      y + data.input_str_dist >= x
    ) {
      data.non_zero_indices[data.non_zero_count++] = 2 * i + offset;
      data.probabilities[2 * i + offset]           = data.solutions[i].real();
      sum += data.solutions[i].real();
    }
  }
  
  // Normalizing pass.
  for (i = 0; i < (int)pow((double)data.row_length, 2); ++i) {
    data.probabilities[i] /= sum;
    normal_sum            += data.probabilities[i];
  }
  
  fftw_destroy_plan(plan);
  #ifdef FFTBOR_DEBUG
  printf("\nPartition function: ");
  printf(data.precision_format, data.partition_function);
  printf("\nSum of eligible data.probabilities > 0: ");
  printf(data.precision_format, sum);
  printf("\nSum of normalized data.probabilities: ");
  printf(data.precision_format, normal_sum);
  printf("\n");
  #endif
}

void print_output(FFTBOR2D_PARAMS& parameters, FFTBOR2D_DATA& data) {
  int i;
  int matrix_size = (int)pow((double)data.row_length, 2);
  
  if (parameters.format == 'B') {
    printf("%s\n%s\n%s\n", parameters.sequence, parameters.structure_1, parameters.structure_2);
    printf("%d,%d\n", data.input_str_dist, data.row_length);
    printf("k\tl\tp(Z_{k,l}/Z)\t-RTln(Z_{k,l})\n");
  }
  
  if (parameters.format == 'M') {
    for (i = 0; i < matrix_size; ++i) {
      if (i && !(i % data.row_length)) {
        printf("\n");
      }
      
      printf(data.precision_format, data.probabilities[i]);
      printf("\t");
    }
    
    printf("\n");
  } else {
    for (i = 0; i < matrix_size; ++i) {
      if (data.probabilities[i] > 0) {
        printf("%d\t%d\t", i / data.row_length, i % data.row_length);
        printf(data.precision_format, data.probabilities[i]);
        printf("\t");
        printf(data.precision_format, -(data.RT / 100) * log(data.probabilities[i]) - (data.RT / 100) * log(data.partition_function));
        printf("\n");
      }
    }
  }
}

inline int j_paired_in(int i, int j, int* base_pairs) {
  return base_pairs[j] >= i && base_pairs[j] < j ? 1 : 0;
}

void populate_matrices(dcomplex* roots_of_unity, int num_roots) {
  int i;
  #pragma omp parallel for default(shared)
  
  for (i = 0; i < num_roots; ++i) {
    roots_of_unity[i] = dcomplex(cos(-2 * M_PI * i / num_roots), sin(-2 * M_PI * i / num_roots));
  }
}

inline void flush_matrices(dcomplex** Z, dcomplex** ZB, dcomplex** ZM, dcomplex** ZM1, int sequence_length) {
  int i, j;
  
  for (i = 0; i <= sequence_length; ++i) {
    for (j = 0; j <= sequence_length; ++j) {
      if (i > 0 && j > 0 && abs(j - i) <= MIN_PAIR_DIST) {
        Z[i][j] = ONE_C;
      } else {
        Z[i][j] = ZERO_C;
      }
      
      ZB[i][j]  = ZERO_C;
      ZM[i][j]  = ZERO_C;
      ZM1[i][j] = ZERO_C;
    }
  }
}

int* get_bp_list(char* sec_str) {
  /* Returns list L of ordered pairs (i,j) where i<j and
   * positions i,j occupied by balancing parentheses
   * For linear time efficiency, use stack
   * Assume that secStr is string consisting of '(',')' and '.'
   * Values -2,-1 returned mean NOT well balanced
   * -2 means too many ( with respect to )
   * -1 means too many ) with respect to (
   * If 1,-1 not returned, then return (possibly empty) list */
  int len = strlen(sec_str);
  int* s = (int*) calloc(len / 2, sizeof(int)); //empty stack
  int* l = (int*) calloc(2 * len * (len - 1) / 2 + 1, sizeof(int)); /* initially empty
                   * list of base
                   * pairs */
  int j, k = 0;
  char ch;
  /* First position holds the number of base pairs */
  l[0] = 0;
  
  for (j = 1; j <= len; j++) {
    l[j] = -1;
  }
  
  for (j = 1; j <= len; j++) {
    ch = sec_str[j - 1];
    
    if (ch == '(') {
      s[k++] = j;
    } else if (ch == ')') {
      if (k == 0) {
        /* There is something wrong with the structure. */
        l[0] = -1;
        return l;
      } else {
        l[s[--k]] = j;
        l[j] = s[k];
        l[0]++;
      }
    }
  }
  
  if (k != 0) {
    /* There is something wrong with the structure. */
    l[0] = -2;
  }
  
  free(s);
  return l;
}

/* Number of base pairs in the region i to j in bpList */
int num_bp(int i, int j, int* bp_list) {
  int n = 0;
  int k;
  
  for (k = i; k <= j; k++)
    if (k < bp_list[k] && bp_list[k] <= j) {
      n++;
    }
    
  return n;
}

short rna_int_code(char a) {
  /* Return the integer corresponding to the given base
     @ = 0  A = 1  C = 2  G = 3  U = 4 */
  if (a == 'A') {
    return 1;
  } else if (a == 'C') {
    return 2;
  } else if (a == 'G') {
    return 3;
  } else if (a == 'U') {
    return 4;
  } else {
    return 0;
  }
}

void initialize_can_base_pair_matrix(int** can_base_pair) {
  // A = 1, C = 2, G = 3, U = 4
  can_base_pair[1][4] = 5;
  can_base_pair[4][1] = 6;
  can_base_pair[2][3] = 1;
  can_base_pair[3][2] = 2;
  can_base_pair[3][4] = 3;
  can_base_pair[4][3] = 4;
}

void translate_to_int_sequence(char* a, short* int_sequence) {
  int i;
  int_sequence[0] = strlen(a);
  
  for (i = 0; i < int_sequence[0]; ++i) {
    int_sequence[i + 1] = rna_int_code(a[i]);
  }
}

void initialize_base_pair_count_matrix(int** num_base_pairs, int* bp_list, int n) {
  int d, i, j;
  
  for (i = 1; i <= n; ++i) {
    for (j = 1; j <= n; ++j) {
      num_base_pairs[i][j] = 0;
    }
  }
  
  for (d = MIN_PAIR_DIST + 1; d < n; d++) {
    for (i = 1; i <= n - d; ++i) {
      j = i + d;
      num_base_pairs[i][j] = num_bp(i, j, bp_list);
    }
  }
}

// INLINE  PRIVATE int E_Hairpin(int size, int type, int si1, int sj1, const char *string, paramT *P){
inline double hairpin_loop_energy(FFTBOR2D_DATA& data, int i, int j, int type, short si1, short sj1, char* string) {
  double energy;
  int size = j - i - 1;
  energy = (size <= 30) ? data.vienna_params->hairpin[size] : data.vienna_params->hairpin[30] + (int)(data.vienna_params->lxc * log((size) / 30.));
  
  if (data.vienna_params->model_details.special_hp) {
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7] = {0}, *ts;
      strncpy(tl, string, 6);
      
      if ((ts = strstr(data.vienna_params->Tetraloops, tl))) {
        return (data.vienna_params->Tetraloop_E[(ts - data.vienna_params->Tetraloops) / 7]);
      }
    } else if (size == 6) {
      char tl[9] = {0}, *ts;
      strncpy(tl, string, 8);
      
      if ((ts = strstr(data.vienna_params->Hexaloops, tl))) {
        return (energy = data.vienna_params->Hexaloop_E[(ts - data.vienna_params->Hexaloops) / 9]);
      }
    } else if (size == 3) {
      char tl[6] = {0, 0, 0, 0, 0, 0}, *ts;
      strncpy(tl, string, 5);
      
      if ((ts = strstr(data.vienna_params->Triloops, tl))) {
        return (data.vienna_params->Triloop_E[(ts - data.vienna_params->Triloops) / 6]);
      }
      
      return (energy + (type > 2 ? data.vienna_params->TerminalAU : 0));
    }
  }
  
  energy += data.vienna_params->mismatchH[type][si1][sj1];
  return energy;
}

// INLINE  PRIVATE int E_IntLoop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P){
inline double interior_loop_energy(FFTBOR2D_DATA& data, int i, int j, int k, int l, int type, int type_2, short si1, short sq1, short sj1, short sp1) {
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns;
  double energy;
  int n1 = k - i - 1;
  int n2 = j - l - 1;
  energy = INF;
  
  if (n1 > n2) {
    nl = n1;
    ns = n2;
  } else {
    nl = n2;
    ns = n1;
  }
  
  if (nl == 0) {
    return data.vienna_params->stack[type][type_2];  /* stack */
  }
  
  if (ns == 0) {                    /* bulge */
    energy = (nl <= MAXLOOP) ? data.vienna_params->bulge[nl] :
             (data.vienna_params->bulge[30] + (int)(data.vienna_params->lxc * log(nl / 30.)));
             
    if (nl == 1) {
      energy += data.vienna_params->stack[type][type_2];
    } else {
      if (type > 2) {
        energy += data.vienna_params->TerminalAU;
      }
      
      if (type_2 > 2) {
        energy += data.vienna_params->TerminalAU;
      }
    }
    
    return energy;
  } else {                          /* interior loop */
    if (ns == 1) {
      if (nl == 1) {                /* 1x1 loop */
        return data.vienna_params->int11[type][type_2][si1][sj1];
      }
      
      if (nl == 2) {                /* 2x1 loop */
        if (n1 == 1) {
          energy = data.vienna_params->int21[type][type_2][si1][sq1][sj1];
        } else {
          energy = data.vienna_params->int21[type_2][type][sq1][si1][sp1];
        }
        
        return energy;
      } else { /* 1xn loop */
        energy = (nl + 1 <= MAXLOOP) ? (data.vienna_params->internal_loop[nl + 1]) : (data.vienna_params->internal_loop[30] + (int)(data.vienna_params->lxc * log((nl + 1) / 30.)));
        energy += MIN2(MAX_NINIO, (nl - ns) * data.vienna_params->ninio[2]);
        energy += data.vienna_params->mismatch1nI[type][si1][sj1] + data.vienna_params->mismatch1nI[type_2][sq1][sp1];
        return energy;
      }
    } else if (ns == 2) {
      if (nl == 2)      {           /* 2x2 loop */
        return data.vienna_params->int22[type][type_2][si1][sp1][sq1][sj1];
      } else if (nl == 3) {         /* 2x3 loop */
        energy = data.vienna_params->internal_loop[5] + data.vienna_params->ninio[2];
        energy += data.vienna_params->mismatch23I[type][si1][sj1] + data.vienna_params->mismatch23I[type_2][sq1][sp1];
        return energy;
      }
    }
    
    {
      /* generic interior loop (no else here!)*/
      energy = (n1 + n2 <= MAXLOOP) ? (data.vienna_params->internal_loop[n1 + n2]) : (data.vienna_params->internal_loop[30] + (int)(data.vienna_params->lxc * log((n1 + n2) / 30.)));
      energy += MIN2(MAX_NINIO, (nl - ns) * data.vienna_params->ninio[2]);
      energy += data.vienna_params->mismatchI[type][si1][sj1] + data.vienna_params->mismatchI[type_2][sq1][sp1];
    }
  }
  
  return energy;
}
