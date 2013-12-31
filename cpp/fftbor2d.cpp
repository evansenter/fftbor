#include <stdio.h>
#include <sys/time.h>
#include "params.h"
#include "initializers.h"
#include "functions.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define TIMING(start, stop, task) printf("[benchmarking] %8.2f\ttime in ms for %s\n", (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0), task);

int main(int argc, char** argv) {
  struct timeval full_start, full_stop, start, stop;
  FFTBOR2D_PARAMS parameters;
  parameters = parse_fftbor2d_args(argc, argv);
  
  if (parameters.benchmark) {
    gettimeofday(&full_start, NULL);
    gettimeofday(&start, NULL);
  }
  
  FFTBOR2D_DATA data;
  FFTBOR2D_THREADED_DATA* threaded_data;
  data          = init_fftbor2d_data(parameters);
  threaded_data = init_fftbor2d_threaded_data(parameters, data);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "initialization")
    gettimeofday(&start, NULL);
  }
  
  precalculate_energies(data);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "pre-calculate energies")
    gettimeofday(&start, NULL);
  }
  
  evaluate_recursions_in_parallel(parameters, data, threaded_data);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "explicit evaluation of Z")
    gettimeofday(&start, NULL);
  }
  
  // Convert point-value solutions to coefficient form w/ inverse DFT.
  populate_remaining_roots(data);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "inferred evaluation of Z")
    gettimeofday(&start, NULL);
  }
  
  solve_system(parameters, data);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "FFT")
    gettimeofday(&start, NULL);
  }
  
  if (parameters.verbose) {
    print_fftbor2d_data(data);
  }
  
  #ifndef SILENCE_OUTPUT
  print_output(parameters, data);
  #endif
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "print output / matrix inversion (if -x was provided)")
    gettimeofday(&start, NULL);
  }
  
  free_fftbor2d_threaded_data(threaded_data, parameters.max_threads);
  free_fftbor2d_data(data);
  free_fftbor2d_params(parameters);
  
  if (parameters.benchmark) {
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "free memory")
  }
  
  if (parameters.benchmark) {
    gettimeofday(&full_stop, NULL);
    TIMING(full_start, full_stop, "total")
  }
  
  return 0;
}
