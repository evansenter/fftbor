#ifndef LIBMFPT_HEADER_H
#define LIBMFPT_HEADER_H

typedef struct {
  int start_state;
  int end_state;
  int sequence_length;
  short energy_based;
  short transition_matrix_input;
  short pseudoinverse;
  short single_bp_moves_only;
  short hastings;
  short verbose;
  double additive_epsilon;
  double distributed_epsilon;
} MFPT_PARAMETERS;

#ifdef __cplusplus
  extern "C" {
    MFPT_PARAMETERS init_mfpt_params();
    int mfpt_error_handling(MFPT_PARAMETERS);
    double** convert_energy_grid_to_transition_matrix(int**, int**, double**, unsigned long*, MFPT_PARAMETERS);
    double compute_mfpt(int*, int*, double**, unsigned long, MFPT_PARAMETERS);
  }
#else
  extern MFPT_PARAMETERS init_mfpt_params();
  int mfpt_error_handling(MFPT_PARAMETERS);
  extern double** convert_energy_grid_to_transition_matrix(int**, int**, double**, unsigned long*, MFPT_PARAMETERS);
  extern double compute_mfpt(int*, int*, double**, unsigned long, MFPT_PARAMETERS);
#endif

#endif
