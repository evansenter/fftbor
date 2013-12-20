#ifndef LIBSPECTRAL_HEADER_H
#define LIBSPECTRAL_HEADER_H

typedef struct {
  int verbose;
  char* sequence;
  char* energy_grid_file;
  char* start_structure;
  char* end_structure;
  double temperature;
  double start_time;
  double end_time;
  double step_size;
  int lonely_bp;
  int energy_cap;
  int eigen_only;
  int use_min;
} SPECTRAL_PARAMS;

typedef struct {
  double *values;       
  double *vectors;
  double *inverse_vectors;
} EIGENSYSTEM;

#ifdef __cplusplus
  extern "C" {
    SPECTRAL_PARAMS init_spectral_params();
    int spectral_error_handling(SPECTRAL_PARAMS);
    EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double*, int);
    void invert_matrix(EIGENSYSTEM, int);
    double probability_at_time(EIGENSYSTEM, double, int, int, int);
  }
#else
  extern SPECTRAL_PARAMS init_spectral_params();
  extern int spectral_error_handling(SPECTRAL_PARAMS);
  extern EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double*, int);
  extern void invert_matrix(EIGENSYSTEM, int);
  extern double probability_at_time(EIGENSYSTEM, double, int, int, int);
#endif

#endif
