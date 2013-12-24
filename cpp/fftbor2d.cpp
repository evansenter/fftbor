#include <stdlib.h>
#include "params.h"
#include "functions.h"

int main(int argc, char *argv[]) {
  char *sequence;            // sequence points to array a[0], ..., a[n] where a[0] = n, and a[1],...,a[n] are RNA nucleotides and where n <= Nseq - 2
  int **intBP = new int*[2]; // intBP is an array of basepairs, where a base pair is defined by two integers
  
  if (argc == 1) {
    usage();
    exit(1);
  }
  
  read_input(argc, argv, &sequence, intBP);
  neighbors(sequence, intBP);

  free(sequence);
  free(intBP[0]);
  free(intBP[1]);
  delete[] intBP;

  return 0;
}
