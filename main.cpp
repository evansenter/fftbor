#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <iostream>
#include "partition.h"
#include "misc.h"
#include "params.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

int           PF, N, PRECISION, MAXTHREADS, ROW_LENGTH, MATRIX_FORMAT, SIMPLE_OUTPUT;
extern double temperature;
char          *ENERGY;
paramT        *P;

void read_input(int, char **, char **, int **);
void usage();

int main(int argc, char *argv[]) {
  char *a;  /* a points to array a[0], ..., a[n] where a[0] = n, and a[1],...,a[n] are RNA nucleotides and where n <= Nseq - 2 */
  int **bps = new int*[2]; /* bps is an array of basepairs, where a base pair is defined by two integers */
  
  if (argc == 1) {
    usage();
    exit(1);
  }
  
  read_input(argc, argv, &a, bps);
  neighbors(a, bps);

  free(a);
  free(bps[0]);
  free(bps[1]);
  delete []bps;

  return 0;
}

void usage() {
  fprintf(stderr, "FFTbor2D [options] sequence structure_1 structure_2\n\n");
  
  fprintf(stderr, "FFTbor2D [options] filename\n");
  fprintf(stderr, "where filename is a file of the format:\n");
  fprintf(stderr, "\t>comment (optional line)\n");
  fprintf(stderr, "\tsequence\n");
  fprintf(stderr, "\tsecondary structure (1)\n");
  fprintf(stderr, "\tsecondary structure (2)\n\n");
  
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-E\tenergyfile,    the default is rna_turner1999.par in this executable's directory. Must be the name of a file with all energy parameters (in the same format as used in Vienna RNA).\n");
  fprintf(stderr, "-T\ttemperature,   the default is 37 degrees Celsius (unless an energyfile with parameters for a different temperature is used.\n");
  fprintf(stderr, "-P\tprecision,     the default is 8, indicates the precision of the probabilities Z_k / Z to be returned (0-9, 0 disables precision handling).\n");
  fprintf(stderr, "-R\trow length,    the default is 100, takes an integer 0 < R <= 100 that describes the dimensions of the 2D matrix in terms of the percentage sequence length.\n");
  fprintf(stderr, "-M\tmatrix format, the default is disabled, presents output in a matrix format instead of a column format.\n");
  fprintf(stderr, "-S\tsimple output, the default is disabled, presents output in column format, for non-zero entries only with no header output (columns are: k, l, p(Z_{k,l}/Z), -RTln(Z_{k,l})).\n");
  
  exit(1);
}

void read_input(int argc, char *argv[], char **maina, int **bps) {
  FILE *infile;
  char line[MAX_SEQ_LENGTH];
  int i, k;
  char *seq = NULL, *str1 = NULL, *str2 = NULL;
 
  MAXTHREADS    = omp_get_max_threads(); 
  PF            = 0;
  PRECISION     = 4;
  ROW_LENGTH    = 0;
  MATRIX_FORMAT = 0; 
  SIMPLE_OUTPUT = 0; 
  ENERGY        = (char *)"rna_turner1999.par";

  /* Function to retrieve RNA sequence and structure, when
   * either input in command line or in a file, where the first
   * line can be a comment (after a >). */
 
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (strcmp(argv[i], "-T") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%lf", &temperature)) {
          usage();
        }
      } else if (strcmp(argv[i], "-E") == 0) {
        if (i == argc - 1) {
          usage();
        }
        ENERGY = argv[++i];
      } else if (strcmp(argv[i], "-R") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &ROW_LENGTH)) {
          usage();
        } else if (ROW_LENGTH < 0 || ROW_LENGTH > 100) {
          usage();
        }
      } else if (strcmp(argv[i], "-M") == 0) {
        if (i == argc - 1) {
          usage();
        }
        MATRIX_FORMAT = 1;
      } else if (strcmp(argv[i], "-S") == 0) {
        if (i == argc - 1) {
          usage();
        }
        SIMPLE_OUTPUT = 1;
      } else if (strcmp(argv[i], "-P") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &PRECISION)) {
          usage();
        } else if (PRECISION < 0 || PRECISION > 9) {
          usage();
        }
      } else {
        usage();
      }
    } else {
      /* File as input? */
      infile = fopen(argv[i], "r");
      if (infile == NULL) { 
        /* Input is not a file */
	      /* arMAX_SEQ_LENGTH[i] should be the sequence and argv[i+1] should be the structure */
        if (argc <= i + 1) {
          usage();
        }
        N        = strlen(argv[i]);
        (*maina) = (char *)xcalloc(N + 1, sizeof(char));
	      seq      = *maina;
	      str1     = (char *)xcalloc(N + 1, sizeof(char));
	      str2     = (char *)xcalloc(N + 1, sizeof(char));
	      (void)sscanf(argv[i++], "%s", seq);
	      (void)sscanf(argv[i++], "%s", str1);
	      (void)sscanf(argv[i], "%s", str2);
	      N = strlen(seq);
        
        if (!TRIEQUALS(strlen(seq), strlen(str1), strlen(str2))) {
          printf("%s\n%s\n%s\n", seq, str1, str2);
          fprintf(stderr,"Length of RNA sequence and structures must be equal\n");
          exit(1);
        }
        
        /* Convert RNA sequence to uppercase and make sure there are no Ts in the sequence (replace by U) */
        for (k = 0; k < N; k++) {
          seq[k] = toupper(seq[k]);
          if (seq[k] == 'T') {
            seq[k] = 'U';
          }
        }
        seq[N]  = '\0';
        str1[N] = '\0';
        str2[N] = '\0';
      }
      else { 
        /* Input is a file */
        if (fgets(line, sizeof(line), infile) == NULL) {
          fprintf(stderr, "There was an error reading the file\n");
          exit(1);
        }

        while ((*line == '*') || (*line == '\0') || (*line == '>')) {
          if (fgets(line, sizeof(line), infile) == NULL) {
            break;
          }
        } 
        
        if ((line == NULL)) {
          usage();
        }
        
        N        = strlen(line);
        (*maina) = (char *)xcalloc(N + 1, sizeof(char));
        seq      = *maina;
        (void)sscanf(line, "%s", seq);
        
        for (k = 0; k < N; k++) {
          seq[k] = toupper(seq[k]);
          if (seq[k] == 'T') {
            seq[k] = 'U';
          }
        }
	
        if (fgets(line, sizeof(line), infile) == NULL) {
          fprintf(stderr,"There was an error reading the file\n");
          exit(1);
        }

        str1 = (char *)xcalloc(N + 1, sizeof(char));
        (void)sscanf(line, "%s", str1);
        
        if (fgets(line, sizeof(line), infile) == NULL) {
          fprintf(stderr,"There was an error reading the file\n");
          exit(1);
        }
        
        str2 = (char *)xcalloc(N + 1, sizeof(char));
        (void)sscanf(line, "%s", str2);

        if (!TRIEQUALS(strlen(seq), strlen(str1), strlen(str2))) {
          printf("%s\n%s\n%s\n", seq, str1, str2);
          fprintf(stderr, "Length of RNA sequence and structure must be equal\n");
          exit(1);
        }
        
        if (N < (int)strlen(seq)) {
          fprintf(stderr,"Length of RNA exceeds array size %d\n",N);
          exit(1);
        } 
        
        fclose(infile);
      }
    }
  }
  
  /* Post-sequence / structure validations */
  if (seq == NULL || str1 == NULL || str2 == NULL) {
    usage();
  }
  
  if (!SIMPLE_OUTPUT) {
    /* Print sequence length, sequence and starting structure */
    printf("%s\n%s\n%s\n", seq, str1, str2);
  }

    bps[0] = getBasePairList(str1);
    bps[1] = getBasePairList(str2);
} 
