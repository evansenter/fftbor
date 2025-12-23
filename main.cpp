#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <cctype>
#include <new>  // for std::bad_alloc
#include <stdexcept>
#include "delta.h"
#include "misc.h"
#include "params.h"
#include "memory_types.h"

// Global variables - defined in globals.cpp
extern int    PF, N, PRECISION, WINDOW_SIZE, MIN_WINDOW_SIZE;
extern double temperature;
extern char   *ENERGY;
extern fftbor::ParamPtr P;

void read_input(int argc, char *argv[], fftbor::CharSequencePtr& seq, fftbor::BasePairListPtr& bps);
[[noreturn]] void usage();

int main(int argc, char *argv[]) {
  try {
    fftbor::CharSequencePtr seq;  // RNA sequence
    fftbor::BasePairListPtr bps;  // Base pair list

    if (argc == 1) {
      usage();
    }

    read_input(argc, argv, seq, bps);
    neighbours(seq.get(), bps.get());

    // Smart pointers automatically clean up - no manual free() needed
    return 0;
  } catch (const std::runtime_error& e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return 1;
  } catch (const std::bad_alloc& e) {
    fprintf(stderr, "Error: Memory allocation failed: %s\n", e.what());
    return 1;
  }
}

[[noreturn]] void usage() {
  fprintf(stderr, "FFTbor [options] sequence structure [options]\n\n");

  fprintf(stderr, "FFTbor [options] filename [options]\n");
  fprintf(stderr, "where filename is a file of the format:\n");
  fprintf(stderr, "\t>comment (optional line)\n");
  fprintf(stderr, "\tsequence\n");
  fprintf(stderr, "\tsecondary structure\n\n");

  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-E\tenergyfile,  the default is rna_turner2004.par in this executable's directory. Must be the name of a file with all energy parameters (in the same format as used in Vienna RNA).\n");
  fprintf(stderr, "-T\ttemperature, the default is 37 degrees Celsius (unless an energyfile with parameters for a different temperature is used.\n");
  fprintf(stderr, "-P\tprecision,   the default is 4, indicates the precision of the probabilities Z_k / Z to be returned (0-9, 0 disables precision handling).\n");

  throw std::runtime_error("Invalid usage");
}

void read_input(int argc, char *argv[], fftbor::CharSequencePtr& main_seq, fftbor::BasePairListPtr& bps) {
  FILE *infile;
  char line[MAX_SEQ_LENGTH];
  int i, k;
  char *seq = NULL;
  fftbor::CharSequencePtr str_ptr;  // Structure string (local, auto-freed)

  PF              = 0;
  PRECISION       = 4;
  WINDOW_SIZE     = 0;
  MIN_WINDOW_SIZE = 0;
  ENERGY          = (char *)"rna_turner2004.par";

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
	      /* argv[i] should be the sequence and argv[i+1] should be the structure */
        if (argc <= i + 1) {
          usage();
        }
        N = strlen(argv[i]);
        main_seq.reset((char *)xcalloc(N + 1, sizeof(char)));
        seq = main_seq.get();
        str_ptr.reset((char *)xcalloc(N + 1, sizeof(char)));
        (void)sscanf(argv[i++], "%s", seq);
        (void)sscanf(argv[i], "%s", str_ptr.get());
        N = strlen(seq);

        if (strlen(seq) != strlen(str_ptr.get())) {
          throw std::runtime_error("Length of RNA sequence and structure must be equal");
        }

        /* Convert RNA sequence to uppercase and make sure there are no Ts in the sequence (replace by U) */
        for (k = 0; k < N; k++) {
          seq[k] = toupper(seq[k]);
          if (seq[k] == 'T') {
            seq[k] = 'U';
          }
        }
        str_ptr.get()[N] = '\0';
        seq[N] = '\0';
      }
      else { 
        /* Input is a file */
        if (fgets(line,sizeof(line),infile) == NULL) {
          throw std::runtime_error("Error reading input file");
        }

        while ((*line == '*') || (*line == '\0') || (*line == '>')) {
          if (fgets(line, sizeof(line), infile) == NULL) {
            break;
          }
        }

        // Check if line is empty (fgets exhausted or only whitespace/comments)
        if (*line == '\0' || *line == '*' || *line == '>') {
          fprintf(stderr, "Error: No sequence found in input file\n");
          usage();
        }
        
        N = strlen(line);
        main_seq.reset((char *)xcalloc(N+1,sizeof(char)));
        seq = main_seq.get();
        (void)sscanf(line,"%s",seq);
        
        for (k = 0; k < N; k++) {
          seq[k] = toupper(seq[k]);
          if (seq[k] == 'T') {
            seq[k] = 'U';
          }
        }
	
        if (fgets(line, sizeof(line), infile) == NULL) {
          throw std::runtime_error("Error reading structure from input file");
        }

        str_ptr.reset((char *)xcalloc(N + 1, sizeof(char)));
        (void)sscanf(line, "%s", str_ptr.get());

        if (strlen(seq) != strlen(str_ptr.get())) {
          throw std::runtime_error("Length of RNA sequence and structure must be equal");
        }

        if (N < (int)strlen(seq)) {
          throw std::runtime_error("Length of RNA exceeds array size");
        } 
        
        fclose(infile);
      }
    }
  }
  
  /* Post-sequence / structure validations */
  if (seq == NULL || !str_ptr) {
    usage();
  }
  
  N = (int)strlen(seq);
  
  if (!WINDOW_SIZE) {
    WINDOW_SIZE = N;
  }
  
  if (!MIN_WINDOW_SIZE) {
    MIN_WINDOW_SIZE = WINDOW_SIZE;
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
  printf("%d %s %s\n", N, seq, str_ptr.get());

  bps = get_base_pair_list(str_ptr.get());

  /* Check for unbalanced parentheses in structure */
  if (bps[0] < 0) {
    if (bps[0] == -1) {
      throw std::runtime_error("Unbalanced structure - too many ')' parentheses");
    } else if (bps[0] == -2) {
      throw std::runtime_error("Unbalanced structure - too many '(' parentheses");
    } else {
      throw std::runtime_error("Invalid structure notation");
    }
  }

  // str_ptr automatically freed when function returns
} 
