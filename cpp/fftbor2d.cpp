#include <stdlib.h>
#include "params.h"
#include "functions.h"

int main(int argc, char** argv) {
  FFTBOR2D_PARAMS parameters;

  if (argc == 1) {
    fftbor2d_usage();
    exit(0);
  }

  parameters = parse_fftbor2d_args(argc, argv);
  neighbors(parameters);
  
  return 0;
}
