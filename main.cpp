#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "delta.h"
#include "misc.h"

int           PF, N, PRECISION, WINDOW_SIZE, MIN_WINDOW_SIZE;
extern double temperature;
char          *ENERGY;
paramT        *P;

void read_input(int, char **, char **, int**);
void usage();

int main(int argc, char *argv[]) {
  char *a;  /* a points to array a[0],...,a[n] where a[0]=n, and a[1],...,a[n] are RNA nucleotides and where n <= Nseq - 2 */
  int *bps; /* bps is an array of basepairs, where a base pair is defined by two integers */
  
  if (argc == 1) {
    usage();
    exit(1);
  }
  
  read_input(argc, argv, &a, &bps);
  neighbours(a, bps);

  if (PF == 1) {
    printf("Total Z\n");
    pf(a);
  }

  free(a);
  free(bps);

  return 0;
}

void usage() {
  fprintf(stderr, "FFTbor [options] sequence structure [options]\n\n");
  
  fprintf(stderr, "FFTbor [options] filename [options]\n");
  fprintf(stderr, "where filename is a file of the format:\n");
  fprintf(stderr, "\t>comment (optional line)\n");
  fprintf(stderr, "\tsequence\n");
  fprintf(stderr, "\tsecondary structure\n\n");
  
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-pf\tcompute the total partition function Z.\n");
  fprintf(stderr, "-E\tenergyfile,  the default is energy.par in this executable's directory. Must be the name of a file with all energy parameters (in the same format as used in Vienna RNA).\n");
  fprintf(stderr, "-T\ttemperature, the default is 37 degrees Celsius (unless an energyfile with parameters for a different temperature is used.\n");
  fprintf(stderr, "-P\tprecision,   the default is 4, indicates the precision of the probabilities Z_k / Z to be returned (0-9, 0 disables precision handling).\n");
  
  exit(1);
}

void read_input(int argc,char *argv[], char **maina, int **bps) {
  FILE *infile;
  char line[MAX_SEQ_LENGTH];
  int i, k;
  char *seq = NULL, *str = NULL;
  
  PF              = 0;
  PRECISION       = 4;
  WINDOW_SIZE     = 0;
  MIN_WINDOW_SIZE = 0;
  ENERGY          = (char *)"energy.par";

  /* Function to retrieve RNA sequence and structure, when
   * either input in command line or in a file, where the first
   * line can be a comment (after a >). */
 
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (strcmp(argv[i], "-pf") == 0) {
        PF = 1;
      } else if (strcmp(argv[i], "-T") == 0) {
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
      } else if (strcmp(argv[i], "-P") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &PRECISION)) {
          usage();
        } else if (PRECISION < 0 || PRECISION > 9) {
          usage();
        }
      } else if (strcmp(argv[i], "-W") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &WINDOW_SIZE)) {
          usage();
        } else if (WINDOW_SIZE < 0) {
          usage();
        }
      } else if (strcmp(argv[i], "-M") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &MIN_WINDOW_SIZE)) {
          usage();
        } else if (MIN_WINDOW_SIZE <= 0) {
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
	      str      = (char *)xcalloc(N + 1, sizeof(char));
	      (void)sscanf(argv[i++], "%s", seq);
	      (void)sscanf(argv[i], "%s", str);
	      N = strlen(seq);
        
        if (strlen(seq) != strlen(str)) {
          fprintf(stderr,"Length of RNA sequence and structure must be equal\n");
          exit(1);
        }
        
        /* Convert RNA sequence to uppercase and make sure there are no Ts in the sequence (replace by U) */
        for (k = 0; k < N; k++) {
          seq[k] = toupper(seq[k]);
          if (seq[k] == 'T') {
            seq[k] = 'U';
          }
        }
        str[N] = '\0';
        seq[N] = '\0';
      }
      else { 
        /* Input is a file */
        if (fgets(line,sizeof(line),infile) == NULL) {
          fprintf(stderr,"There was an error reading the file\n");
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
        (*maina) = (char *)xcalloc(N+1,sizeof(char));
        seq      = *maina;
        (void)sscanf(line,"%s",seq);
        
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

        str = (char *)xcalloc(N + 1, sizeof(char));
        (void)sscanf(line, "%s", str);

        if (strlen(seq) != strlen(str)) {
          printf("%s\n%s\n", seq, str);
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
  if (seq == NULL || str == NULL) {
    usage();
  }
  
  if (WINDOW_SIZE > N) {
    printf("Error: the window size provided (%d) can't be longer than the input sequence length (%d).\n\n", WINDOW_SIZE, N);
    usage();
  }
  
  if (MIN_WINDOW_SIZE > WINDOW_SIZE) {
    printf("Error: the minimum window size provided (%d) can't be larger than the window size provided (%d).\n\n", MIN_WINDOW_SIZE, WINDOW_SIZE);
    usage();
  }
  
  if (MIN_WINDOW_SIZE > N) {
    printf("Error: the minimum window size provided (%d) can't be larger than the input sequence length (%d).\n\n", MIN_WINDOW_SIZE, N);
    usage();
  }
  
  /* Print sequence length, sequence and starting structure */
  printf("%d %s %s\n", (int)strlen(seq), seq, str);

  *bps = getBasePairList(str);
  free(str);
} 
