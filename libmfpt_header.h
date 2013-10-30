#ifndef LIBMFPT_HEADER_H
#define LIBMFPT_HEADER_H

typedef struct GlobalParameters {
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
} GlobalParameters;

#ifdef __cplusplus
  extern "C" {
    GlobalParameters init_params();
    int error_handling(GlobalParameters);
    double** convert_energy_grid_to_transition_matrix(int**, int**, double**, unsigned long*, GlobalParameters);
    double compute_mfpt(int*, int*, double**, unsigned long, GlobalParameters);
  }
#else
  extern GlobalParameters init_params();
  int error_handling(GlobalParameters);
  extern double** convert_energy_grid_to_transition_matrix(int**, int**, double**, unsigned long*, GlobalParameters);
  extern double compute_mfpt(int*, int*, double**, unsigned long, GlobalParameters);
#endif

#endif
