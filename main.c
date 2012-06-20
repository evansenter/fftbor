/* main.c
 * Compute number of k-neighbours, partition function contribution
 * from k-neighbours or the minimal free energy of the k-neighbours.
 * Input: Either a file consisting of five lines;
 * line 1 should start with a >, after that it can contain anything,
 * line 2 sequence
 * line 3 structure
 * line 4 delta (optional)
 * Or RNAbor [options] sequence structure [options]
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
//#include <sys/time.h>
#include "delta.h"
#include "misc.h"

//#define MAXD

/* Function that reformates the input */
void read_input(int ,char **, char **, int**);
void usage();

int main(int argc, char *argv[]) {

  // variable declarations
  char *a; /*a points to array a[0],...,a[n]
	    *where a[0]=n, and a[1],...,a[n] are RNA nucleotides
	    *and where n <= Nseq-2 */
  int *bps; /*bps is an array of basepairs, where a base pair is
	     *defined by two integers*/

  /*---- abort if no command line parameters given ----*/
  if ( argc == 1) {
    usage();
    exit(1);
    }
 

  /* Read input */
  read_input(argc,argv,&a,&bps);

  if (NUMBER)
#ifdef COMPUTEMFE
    printf("k\tZk\tNk\tEk\tstr\n");
#else
    printf("k\tZk\tEk\tstr\n");
#endif
  else
    printf("k\tZk\tEk\tstr\n");
  neighbours(a, bps);

  if (PF==1) {
    printf("Total Z\n");
    pf(a);
  }

  free(a);
  free(bps);

  return 0;
}

void usage() {
  fprintf(stderr,"Usage:\nRNAbor [options] filename [options]\n");
  fprintf(stderr,"where filename has the format\n");
  fprintf(stderr,"\t>comment (optional line)\n");
  fprintf(stderr,"\tsequence\n");
  fprintf(stderr,"\tsecondary structure\n");
  fprintf(stderr,"\tdelta (optional)\n");
  fprintf(stderr,"or\n");
  fprintf(stderr,"RNAbor [options] sequence structure [options]\n\n");
  fprintf(stderr,"the options are the following;\n");
  fprintf(stderr,"-, makes RNAbor to read sequence and structure from stdin\n");
  fprintf(stderr,"-d delta, where delta is the maximum distance from the input \n");
  fprintf(stderr,"          structure to compute delta-neighbours for.\n");
  fprintf(stderr,"          The default value for delta is 5000, i.e. the default is\n");
  fprintf(stderr,"          to limit the maximum delta with the option -s, see below.\n");
  fprintf(stderr,"-s n, will stop the computation after Zk == 0.0 n times in a row\n");
  fprintf(stderr,"      The default value is n=5\n");
  fprintf(stderr,"-pf, compute the total partition function Z = sum_k(Zk)\n");
  fprintf(stderr,"     This option can be useful if -d is used to stop the computation\n");
  fprintf(stderr,"     before all neighbours are counted, since the full\n");
  fprintf(stderr,"     partition function is necessary for computing probabilities.\n");
  fprintf(stderr,"-nodangle, do not include dangling ends in the energy\n");
  fprintf(stderr,"-N 1, will compute the number of k-neighbours (default)\n");
  fprintf(stderr,"-N 0, will only compute the partition function and structure\n");
  fprintf(stderr,"      of the k-neighbours, not the number\n");
  fprintf(stderr,"-E energyfile, where energyfile is the name of a file\n");
  fprintf(stderr,"                with all energy parameters (in the same\n");
  fprintf(stderr,"                format as used in Vienna RNA).\n");
  fprintf(stderr,"-T temp, set the temperature,\n");
  fprintf(stderr,"         the default is 37 degrees Celsius (unless an\n");
  fprintf(stderr,"         energyfile with parameters for a different\n");
  fprintf(stderr,"         temperature is used.\n");
  exit(1);
}

void read_input(int argc,char *argv[],char **maina, int **bps){
  FILE * infile;
  char line[NN];
  int i,k,r;
  char *seq=NULL, *str=NULL;
  
  PF = 0;
  DELTA = 5000;
  DANGLE = 1;
  NUMBER = 1;
  STRUCTURE = 0;
  PARTITION = 1;
  STOP = 5;
  ENERGY = (char *)"energy.par";

  /* Function to retrieve RNA sequence and structure and delta, 
   * either input in command line or in a file, where the first
   * line can be a comment (after a >). */
 
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') {
      if (argv[i][1]=='d') {
	/* delta value */
	if (argv[i][2]!='\0') {
	  r=sscanf(argv[i]+2, "%d", &DELTA);
	  if (!r) usage();
	}
	else if (i==argc-1) usage();
	else if (!sscanf(argv[++i],"%u", &DELTA)) usage();
      }
      else if (strcmp(argv[i],"-pf")==0)
	PF = 1;
      else if (strcmp(argv[i],"-nodangle")==0)
	DANGLE = 0;
      else if (strcmp(argv[i],"-s")==0) {
	if (i==argc-1) usage();
	else if (!sscanf(argv[++i],"%u", &STOP)) usage();
      }
      else if (strcmp(argv[i],"-T")==0) {
	if (i==argc-1) usage();
	else if (!sscanf(argv[++i],"%lf", &temperature)) usage();
      }
      else if (strcmp(argv[i],"-N")==0) {
	if (i==argc-1) usage();
	else if (!sscanf(argv[++i],"%u", &NUMBER)) usage();
      }
      else if (strcmp(argv[i],"-E")==0) {
	if (i==argc-1) usage();
	ENERGY = argv[++i];
      }
      else if (strcmp(argv[i],"-")==0) {
	/* Read from stdin */
	if (fgets(line,sizeof(line),stdin) == NULL) usage();
	while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	  if (fgets(line,sizeof(line),stdin) == NULL) break;
	} 
	if ((line ==NULL) ) usage();
	N = strlen(line);
	(*maina) = (char *) xcalloc(N+1,sizeof(char));
	seq = *maina;
	(void) sscanf(line,"%s",seq);
	for (k = 0;k<N;k++) {
	  seq[k] = toupper(seq[k]);
	  if (seq[k] == 'T') seq[k] = 'U';
	}
	if (fgets(line,sizeof(line),stdin) == NULL) {
	  fprintf(stderr,"There was an error reading from stdin\n");
	  exit(1);
	}
	str = (char *) xcalloc(N+1,sizeof(char));
	(void) sscanf(line,"%s",str);
	if (fgets(line,sizeof(line),stdin)) {
	  /* Delta is given in the input */
	  DELTA = atoi(line);
	}
	if (strlen(seq)!=strlen(str)) {
	  printf("%s\n%s\n",seq,str);
	  fprintf(stderr,"Length of RNA sequence and structure must be equal\n");
	  exit(1);
	}
	
      }
      else
	usage();
    }
    else {
      /* File as input ?*/
      infile = fopen(argv[i],"r"); // Bug fix 27/1
      if (infile==NULL) { /* Input is not a file */
	/* argv[i] should be the sequence and argv[i+1] should be the
	 * structure*/
	if (argc <= i+1)
	  usage();
	N = strlen(argv[i]);
	(*maina) = (char *) xcalloc(N+1,sizeof(char));
	seq = *maina;
	str =  (char *) xcalloc(N+1,sizeof(char));
	(void) sscanf(argv[i++],"%s",seq);
	(void) sscanf(argv[i],"%s",str);
	N  = strlen(seq);
	if (strlen(seq) != strlen(str)) {
	  fprintf(stderr,"Length of RNA sequence and structure must be equal\n");
	  exit(1);
	}
	/* Convert RNA sequence to uppercase and make sure there are
	 *  no Ts in the sequence (replace by U) 
	 */
	for (k=0;k<N;k++) {
	  seq[k] = toupper(seq[k]);
	  if (seq[k] == 'T') seq[k] = 'U';
	}
	str[N] = '\0';
	seq[N] = '\0';
      }
      else { /* Input is a file */
	if (fgets(line,sizeof(line),infile) == NULL) {
	  fprintf(stderr,"There was an error reading the file\n");
	  exit(1);
	}
	while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	  if (fgets(line,sizeof(line),infile) == NULL) break;
	} 
	if ((line ==NULL) ) usage();
	N = strlen(line);
	(*maina) = (char *) xcalloc(N+1,sizeof(char));
	seq = *maina;
	(void) sscanf(line,"%s",seq);
	for (k = 0;k<N;k++) {
	  seq[k] = toupper(seq[k]);
	  if (seq[k] == 'T') seq[k] = 'U';
	}
	
	if (fgets(line,sizeof(line),infile) == NULL) {
	  fprintf(stderr,"There was an error reading the file\n");
	  exit(1);
	}
	str = (char *) xcalloc(N+1,sizeof(char));
	(void) sscanf(line,"%s",str);

	if (fgets(line,sizeof(line),infile)) {
	  /* Delta is given in the input */
	  DELTA = atoi(line);
	}

	if (strlen(seq)!=strlen(str)) {
	  printf("%s\n%s\n",seq,str);
	  fprintf(stderr,"Length of RNA sequence and structure must be equal\n");
	  exit(1);
	}
	if (N<(int) strlen(seq)){
	  fprintf(stderr,"Length of RNA exceeds array size %d\n",N);
	  exit(1);
	} 
	fclose(infile);
      }
    }
  }
  if (seq==NULL || str==NULL) {
    usage();
  }

  if (N<DELTA)
    DELTA = N; /* To avoid using to much memory*/
  
  /* Print sequence length, sequence and starting structure and
   * dangle/nodangle */
  printf("%d %s %s ",N,seq,str);
  if (DANGLE)
    printf("dangle\n");
  else
    printf("nodangle\n");

  *bps = getBasePairList(str);
  free(str);
} 
