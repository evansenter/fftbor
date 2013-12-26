#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <cctype>
#include <limits>
#include "params.h"
#include "vienna/data_structures.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define TRIEQUALS(x, y, z) ((x == y) && (y == z)) /* Transitivity (BOOM) */

extern double temperature;

FFTBOR2D_PARAMS init_fftbor2d_params() {
  FFTBOR2D_PARAMS parameters = {
    NULL,                                  // sequence
    NULL,                                  // structure_1
    NULL,                                  // structure_2
    NULL,                                  // energy_file
    0,                                     // seq_length
    (int)ceil(log(pow(10., 8)) / log(2.)), // precision
    0,                                     // max_threads
    'B',                                   // format
    0                                      // verbose
  };
  #ifdef _OPENMP
  parameters.max_threads = omp_get_max_threads();
  #endif
  return parameters;
}

void free_fftbor2d_params(FFTBOR2D_PARAMS parameters) {
  free(parameters.sequence);
  free(parameters.structure_1);
  free(parameters.structure_2);
  free(parameters.energy_file);
}

FFTBOR2D_PARAMS parse_fftbor2d_args(int argc, char** argv) {
  int i;
  FFTBOR2D_PARAMS parameters;

  if (argc < 1) {
    fftbor2d_usage();
  }

  parameters = init_fftbor2d_params();

  for (i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      if (!strcmp(argv[i], "-T")) {
        if (i == argc - 1) {
          fftbor2d_usage();
        } else if (!sscanf(argv[++i], "%lf", &temperature)) {
          fftbor2d_usage();
        }
      } else if (!strcmp(argv[i], "-E")) {
        if (i == argc - 1) {
          fftbor2d_usage();
        }

        parameters.energy_file = argv[++i];
      } else if (!strcmp(argv[i], "-P")) {
        if (i == argc - 1) {
          fftbor2d_usage();
        } else if (!sscanf(argv[++i], "%d", &parameters.precision)) {
          fftbor2d_usage();
        } else if (parameters.precision < 0 || parameters.precision > std::numeric_limits<double>::digits) {
          fftbor2d_usage();
        }
      } else if (!strcmp(argv[i], "-M")) {
        parameters.format = 'M';
      } else if (!strcmp(argv[i], "-S")) {
        parameters.format = 'S';
      } else if (!strcmp(argv[i], "-X")) {
        parameters.format = 'X';
      } else if (!strcmp(argv[i], "-Y")) {
        parameters.format = 'Y';
      } else if (!strcmp(argv[i], "-V")) {
        parameters.verbose = 1;
      } else {
        fftbor2d_usage();
      }
    } else {
      parse_fftbor2d_sequence_data(argc, argv, i, parameters);
    }
  }

  if (parameters.energy_file == NULL) {
    parameters.energy_file = find_energy_file((char*)"rna_turner2004.par");
  }

  if (parameters.verbose) {
    debug_fftbor2d_parameters(parameters);
  }

  if (fftbor2d_error_handling(parameters)) {
    fftbor2d_usage();
  }

  if (parameters.format == 'B') {
    printf("%s\n%s\n%s\n", parameters.sequence, parameters.structure_1, parameters.structure_2);
  }

  return parameters;
}

void parse_fftbor2d_sequence_data(int argc, char** argv, int& i, FFTBOR2D_PARAMS& parameters) {
  int j;
  FILE* file;
  char line[1024];
  file = fopen(argv[i], "r");

  if (file == NULL) {
    /* Input is not a file */
    /* argv[i] should be the sequence and argv[i + 1], argv[i + 2] should be the structures */
    if (argc <= i + 2) {
      fftbor2d_usage();
    }

    parameters.seq_length  = strlen(argv[i]);
    parameters.sequence    = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_1 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_2 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    sscanf(argv[i++], "%s", parameters.sequence);
    sscanf(argv[i++], "%s", parameters.structure_1);
    sscanf(argv[i],   "%s", parameters.structure_2);
  } else {
    /* Input is a file */
    if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "There was an error reading the file\n");
      exit(0);
    }

    while (*line == '*' || *line == '\0' || *line == '>') {
      if (fgets(line, sizeof(line), file) == NULL) {
        break;
      }
    }

    if ((line == NULL)) {
      fftbor2d_usage();
    }

    parameters.seq_length  = strlen(line);
    parameters.sequence    = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_1 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    parameters.structure_2 = (char*)calloc(parameters.seq_length + 1, sizeof(char));
    sscanf(line, "%s", parameters.sequence);

    if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "There was an error reading the file\n");
      exit(0);
    }

    sscanf(line, "%s", parameters.structure_1);

    if (fgets(line, sizeof(line), file) == NULL) {
      fprintf(stderr, "There was an error reading the file\n");
      exit(0);
    }

    sscanf(line, "%s", parameters.structure_2);
    fclose(file);
  }

  parameters.sequence[parameters.seq_length]    = '\0';
  parameters.structure_1[parameters.seq_length] = '\0';
  parameters.structure_2[parameters.seq_length] = '\0';

  /* Convert RNA sequence to uppercase and make sure there are no Ts in the sequence (replace by U). */
  for (j = 0; j < parameters.seq_length; ++j) {
    parameters.sequence[j] = toupper(parameters.sequence[j]);

    if (parameters.sequence[j] == 'T') {
      parameters.sequence[j] = 'U';
    }
  }
}

char* find_energy_file(char* energy_file_name) {
  char* envPath, *tempPath, *splitPath, *energyLocation, *possiblePath;
  tempPath       = getenv("PATH");
  envPath        = (char*)calloc((strlen(tempPath) + 2), sizeof(char));
  energyLocation = (char*)calloc(1, sizeof(char));
  strcpy(envPath, ".:");
  strcat(envPath, tempPath);

  if (envPath != NULL) {
    splitPath = strtok(envPath, ":");

    while (splitPath != NULL) {
      possiblePath = (char*)calloc((strlen(splitPath) + 1 + strlen(energy_file_name)), sizeof(char));
      strcpy(possiblePath, splitPath);
      strcat(possiblePath, "/");
      strcat(possiblePath, energy_file_name);

      if (access(possiblePath, R_OK) != -1) {
        energyLocation = (char*)calloc(strlen(possiblePath), sizeof(char));
        strcpy(energyLocation, possiblePath);
        splitPath = NULL;
      } else {
        splitPath = strtok(NULL, ":");
      }

      free(possiblePath);
    }
  }

  free(envPath);
  return energyLocation;
}

int fftbor2d_error_handling(FFTBOR2D_PARAMS parameters) {
  int error = 0;

  if (!TRIEQUALS(strlen(parameters.sequence), strlen(parameters.structure_1), strlen(parameters.structure_2))) {
    fprintf(
      stderr,
      "(%d) %s\n(%d) %s\n(%d) %s\n",
      (int)strlen(parameters.sequence),
      parameters.sequence,
      (int)strlen(parameters.structure_1),
      parameters.structure_1,
      (int)strlen(parameters.structure_2),
      parameters.structure_2
    );
    fprintf(stderr, "Length of RNA sequence and structures must be equal.\n");
    error++;
  }

  if (error) {
    fprintf(stderr, "\n");
  }

  return error;
}

void debug_fftbor2d_parameters(FFTBOR2D_PARAMS parameters) {
  printf("sequence\t\t%s\n",           parameters.sequence    == NULL ? "*missing*" : parameters.sequence);
  printf("structure_1\t\t%s\n",        parameters.structure_1 == NULL ? "*missing*" : parameters.structure_1);
  printf("structure_2\t\t%s\n",        parameters.structure_2 == NULL ? "*missing*" : parameters.structure_2);
  printf("seq_length\t\t%d\n",         parameters.seq_length);
  printf("precision\t\t%d\n",          parameters.precision);
  printf("max_threads\t\t%d\n",        parameters.max_threads);
  printf("energy_file\t\t%s\n",        parameters.energy_file);
  printf("format\t\t%c\n",             parameters.format);
  printf("\n");
}

void fftbor2d_usage() {
  fprintf(stderr, "FFTbor2D [options] sequence structure_1 structure_2\n\n");
  fprintf(stderr, "FFTbor2D [options] filename\n");
  fprintf(stderr, "where filename is a file of the format:\n");
  fprintf(stderr, "\t>comment (optional line)\n");
  fprintf(stderr, "\tsequence (max length: 1024)\n");
  fprintf(stderr, "\tsecondary structure (1)\n");
  fprintf(stderr, "\tsecondary structure (2)\n\n");
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-E\tenergyfile,    the default is rna_turner2004.par in this current directory. Must be the name of a file with all energy parameters (in the same format as used in Vienna RNA). Energy file lookup first checks the current directory, and then iterates through the PATH shell variable until a matching file is found. If no file is found, the default ViennaRNA parameters are used and a warning is presented to the user. If the -E switch is explicitly provided, that file is used in lieu of searching for the rna_turner2004.par file.\n");
  fprintf(stderr, "-T\ttemperature,   the default is 37 degrees Celsius (unless an energyfile with parameters for a different temperature is used.\n");
  fprintf(stderr, "-P\tprecision,     the default is %d, indicates the precision (base 2) of the probabilities Z_k / Z to be returned (0-%d, 0 disables precision handling).\n", (int)ceil(log(pow(10., 8)) / log(2.)), std::numeric_limits<double>::digits);
  fprintf(stderr, "-M\tmatrix format, the default is disabled, presents output in a matrix format instead of a column format.\n");
  fprintf(stderr, "-S\tsimple output, the default is disabled, presents output in column format, for non-zero entries only with no header output (columns are: k, l, p(Z_{k,l}/Z), -RTln(Z_{k,l})).\n");
  fprintf(stderr, "-V\tverbose, the default is disabled, presents some debug information at runtime.\n");
  fprintf(stderr, "-X\tMFPT, the default is disabled, estimates the mean-first passage time of the input RNA from structure 1 to structure 2.\n");
  fprintf(stderr, "-Y\tspectral decomposition, the default is disabled, estimates the population proportion of the input structures over time.\n\n");
  fprintf(stderr, "Note: the output formatting flags (M, S, X, Y) are mutually exclusive. If more than one is provided, *only* the last flag will be honored.\n");
  exit(0);
}
