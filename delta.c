/* delta.c
 * Compute number of delta neighbours and partition functions
 * for different values of k
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include "delta.h"
#include "misc.h"
#include <fftw3.h>
/* Vienna RNA h-files */
#include "energy_const.h"
#include "energy_par.h"
#include <iostream>
#define ZZ(k,i,j) Z[i][j][k]
#define ZZb(k,i,j) Zb[i][j][k]
#define ZZm(k,i,j) Zm[i][j][k]
#ifdef COMPUTEMFE
#define MFEMFE(k,i,j) Z[j][i-1][k]
#define MFEMFEb(k,i,j) Zb[j][i-1][k]
#define MFEMFEm(k,i,j) Zm[j][i-1][k]
#endif
#define NNB(k,i,j) NB[i][j][k]
#define Linterior 30

// Stuff pulled over from rnaborcpp
#define STRUCTURE_COUNT 1
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define PRECISION 4
#define FFTW_REAL 0
#define FFTW_IMAG 1
#define FFTBOR_DEBUG 0

extern int DELTA;
extern int PF;
extern int DANGLE;
extern int NUMBER;
extern char *ENERGY;
extern int STRUCTURE;
extern int PARTITION;
extern int STOP;
extern int N;
extern double temperature;
extern paramT *P;

/* Some Vienna RNA things */

extern "C" void read_parameter_file(const char fname[]);

void neighbours(char *,int *);
void pf(char *);
void sub_structure(int , int , int *, int *);
int numbp(int , int , int *);
void print_bps(int *);
void printbpsstring(int *);
int bn(char );

void initialize_PN(int[][5]);
void translateseq(char *, int *);
void initialize_NumBP(int **, int *, int );
double hairpinloop(int, int, int, int *, char *);
double interiorloop(int, int, int, int, int, int, int *);
double multiloop_closing(int, int, int, int, int, int, int *);

struct index {
  short delta;
  short i;
  short j;
  short type;
};
void backtrack(struct index, struct index ***, struct index ***, struct index ***, int *);
void backtrackb(struct index, struct index ***, struct index ***, struct index ***, int *);
void backtrackm(struct index, struct index ***, struct index ***, struct index ***, int *);

void neighbours(char *a,int *bps) {
  int i, j, k, l, d, delta;
  int n = strlen(a);
  /* Backtrack */
  double RT = 0.0019872370936902486 * (temperature + 273.15) * 100; // 0.01 * (kcal K)/mol

  double energy;

  int *seq = (int *)xcalloc(n+1,sizeof(int));
  static int PN[5][5];
  static int **NumBP;

  const char *fname = ENERGY;
  N=n;
  read_parameter_file(fname);
  P = scale_parameters();

  translateseq(a,seq);
  initialize_PN(PN);
  NumBP = (int **) xcalloc(n+1,sizeof(int *));
  for (i=1;i<=n;i++)
    NumBP[i] = (int *) xcalloc(n+1,sizeof(int));
  initialize_NumBP(NumBP, bps, n);
  
  // ****************************************************************************
  // FFTbor code starts
  // ****************************************************************************
  // Variable declarations.
  int root, **bpCounts;
  double scalingFactor;
  dcomplex x;
  
  int sequenceLength = strlen(a);
  char *sequence     = new char[sequenceLength + 1];

	sequence[0] = '@';
  strncpy(sequence + 1, a, sequenceLength);
  
  dcomplex **Z            = new dcomplex*[sequenceLength + 1];
  dcomplex **ZB           = new dcomplex*[sequenceLength + 1];
  dcomplex **ZM           = new dcomplex*[sequenceLength + 1];
  dcomplex **rootsOfUnity = new dcomplex*[sequenceLength + 1];
  double    *coefficients = new double[sequenceLength + 1];
  
  // Matrix allocation.
  for (i = 0; i <= sequenceLength; ++i) {
    Z[i]               = new dcomplex[sequenceLength + 1];
    ZB[i]              = new dcomplex[sequenceLength + 1];
    ZM[i]              = new dcomplex[sequenceLength + 1];
    rootsOfUnity[i]    = new dcomplex[2];
    rootsOfUnity[i][0] = dcomplex(cos(2 * M_PI * i / (sequenceLength + 1)), sin(2 * M_PI * i / (sequenceLength + 1)));
  }
	
	for (i = 1; i <= sequenceLength; ++i) {
		for (j = 1; j <= sequenceLength; ++j) {
			printf("%d %d: %d\n", i, j, PN[seq[i]][seq[j]]);
		}
	}
	
  bpCounts = fillBasePairCounts(bps, sequenceLength);
  
  // Start main recursions (root <= round(sequenceLength / 2.0) is an optimization for roots of unity).
  for (root = 0; root <= round(sequenceLength / 2.0); ++root) {
    // Flush the matrices.
    for (i = 0; i <= sequenceLength; ++i) {
      for (j = 0; j <= sequenceLength; ++j) {
        Z[i][j]  = ZERO_C;
        ZB[i][j] = ZERO_C;
        ZM[i][j] = ZERO_C;
      }
    }
    
    for (d = 0; d <= MIN_PAIR_DIST; ++d) {
      for (i = 1; i <= sequenceLength - d; ++i) {
        j = i + d;
        
        Z[i][j] = ONE_C;
        
        if (STRUCTURE_COUNT && i != j) {
          Z[j][i] = ONE_C;
        }
      }
    }
    
    x = rootsOfUnity[root][0];
    
    // ****************************************************************************
    // Main recursions
    // ****************************************************************************
    for (d = MIN_PAIR_DIST + 1; d < sequenceLength; ++d) {
      for (i = 1; i <= sequenceLength - d; ++i) {
        j = i + d;
      
        if (canBasePair(i, j, sequence)) {
          // ****************************************************************************
          // Solve ZB 
          // ****************************************************************************
          // In a hairpin, [i + 1, j - 1] unpaired.
          energy    = hairpinloop(i, j, PN[seq[i]][seq[j]], seq, a);
          delta     = bpCounts[i][j] + jPairedTo(i, j, bps);
          ZB[i][j] += pow(x, delta) * exp(-energy / RT);

          if (STRUCTURE_COUNT) {
            ZB[j][i] += 1;
          }

          // Interior loop / bulge / stack / multiloop.
          for (k = i + 1; k < j - MIN_PAIR_DIST; ++k) {
            for (l = max2(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
              if (canBasePair(k, l, sequence)) {
                 // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1) are all unpaired.
                 energy    = interiorloop(i, j, k, l, PN[seq[i]][seq[j]], PN[seq[l]][seq[k]], seq);
                 delta     = bpCounts[i][j] - bpCounts[k][l] + jPairedTo(i, j, bps);
                 ZB[i][j] += (ZB[k][l] * pow(x, delta) * exp(-energy / RT));

                 if (STRUCTURE_COUNT) {
                   ZB[j][i] += ZB[l][k];
                 }

                if (k > i + MIN_PAIR_DIST + 2) {
                  // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, there is at least one hairpin between (i + 1, k - 1).
                  energy    = multiloop_closing(i, j, k, l, PN[seq[j]][seq[i]], PN[seq[k]][seq[l]], seq);
                  delta     = bpCounts[i][j] - bpCounts[i + 1][k - 1] - bpCounts[k][l] + jPairedTo(i, j, bps);
                  ZB[i][j] += ZM[i + 1][k - 1] * ZB[k][l] * pow(x, delta) * exp(-energy / RT);;

                  if (STRUCTURE_COUNT) {
                    ZB[j][i] += ZM[k - 1][i + 1] * ZB[l][k];
                  }
                }
              }
            }
          }
        }
        
        // ****************************************************************************
        // Solve ZM
        // ****************************************************************************
        energy    = P->MLbase;
        delta     = jPairedIn(i, j, bps);
        ZM[i][j] += ZM[i][j - 1] * pow(x, delta) * exp(-energy / RT);

        if (STRUCTURE_COUNT) {
          ZM[j][i] += ZM[j - 1][i];
        }

        for (k = i; k < j - MIN_PAIR_DIST; ++k) {
          if (canBasePair(k, j, sequence)) {
            // Only one stem.
            energy    = P->MLintern[PN[seq[k]][seq[j]]] + P->MLbase * (k - i);
            delta     = bpCounts[i][j] - bpCounts[k][j];
            ZM[i][j] += ZB[k][j] * pow(x, delta) * exp(-energy / RT);

            if (STRUCTURE_COUNT) {
              ZM[j][i] += ZB[j][k];
            }

            // More than one stem.
            if (k > i + THRESHOLD + 2) {
              energy    = P->MLintern[PN[seq[k]][seq[j]]];
              delta     = bpCounts[i][j] - bpCounts[i][k - 1] - bpCounts[k][j];
              ZM[i][j] += ZM[i][k - 1] * ZB[k][j] * pow(x, delta) * exp(-energy / RT);

              if (STRUCTURE_COUNT) {
                ZM[j][i] += ZM[k - 1][i] * ZB[j][k];
              }
            }
          }
        }
        
        // **************************************************************************
        // Solve Z
        // **************************************************************************
        delta    = jPairedIn(i, j, bps);
        Z[i][j] += Z[i][j - 1] * pow(x, delta);

        if (STRUCTURE_COUNT) {
          Z[j][i] += Z[j - 1][i];
        }

        for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
          // (k, j) is the rightmost base pair in (i, j).
          if (canBasePair(k, j, sequence)) {
            energy = PN[seq[k]][seq[j]] > 2 ? TerminalAU : 0;
            
            if (k == i) {
              delta    = bpCounts[i][j] - bpCounts[k][j];
              Z[i][j] += ZB[k][j] * pow(x, delta) * exp(-energy / RT);

              if (STRUCTURE_COUNT) {
                Z[j][i] += ZB[j][k];
              }
            } else {
              delta    = bpCounts[i][j] - bpCounts[i][k - 1] - bpCounts[k][j];
              Z[i][j] += Z[i][k - 1] * ZB[k][j] * pow(x, delta) * exp(-energy / RT);

              if (STRUCTURE_COUNT) {
                Z[j][i] += Z[k - 1][i] * ZB[j][k];
              }
            }
          }
        }
      }
    }
    
    rootsOfUnity[root][1] = Z[1][sequenceLength];
    
    if (!root) {
      scalingFactor = Z[1][sequenceLength].real();
			
			if (FFTBOR_DEBUG) {
				for (k = 0; k <= sequenceLength; ++k) {
					for (l = 0; l <= sequenceLength; ++l) {
						printf("%+9.2f, ", Z[k][l].real());
					}
					printf("\n");
				}
				printf("\n\n");
			}
    }

		if (FFTBOR_DEBUG) {
			std::cout << "." << std::flush;
		}
  }

  // Optimization leveraging complementarity of roots of unity.
  if (sequenceLength % 2) {
    i = root - 2;
  } else {
    i = root - 1;
  }
  
  for (; root <= sequenceLength && i > 0; --i, ++root) {
    rootsOfUnity[root][1] = dcomplex(rootsOfUnity[i][1].real(), -rootsOfUnity[i][1].imag());
  }

	printf("Number of structures: %.0f\n", Z[sequenceLength][1].real());

  solveSystem(sequenceLength, rootsOfUnity, coefficients, scalingFactor);
  // ****************************************************************************
  // FFTbor code ends
  // ****************************************************************************
  
  free(seq);
}

void backtrack(struct index curr, struct index ***T,
	       struct index ***Tb, struct index ***Tm, int *bps) {
  struct index next = T[curr.i][curr.j][curr.delta];
  int d1,d2;
  if (next.type==1) {
    backtrackb(next,T,Tb,Tm,bps);
  }
  else if (next.type==2) {
    d1 = next.delta & 0xff, d2 = (next.delta>>8);
    next.delta = d1;
    backtrackb(next,T,Tb,Tm,bps);
    next.j=next.i-1;
    next.i=curr.i;
    next.delta = d2;
    backtrack(next,T,Tb,Tm,bps);
  }
  else if (next.type==5) {
    backtrack(next,T,Tb,Tm,bps);
  }
}

void backtrackb(struct index curr, struct index ***T, 
		struct index ***Tb, struct index ***Tm, int *bps) {
  struct index next = Tb[curr.i][curr.j][curr.delta];
  int d1,d2;
  bps[0]++;
  bps[curr.i] = curr.j;
  bps[curr.j] = curr.i;
  if (next.type==1) {
    backtrackb(next,T,Tb,Tm,bps);
  }
  else if (next.type==2) {
    d1 = next.delta & 0xff, d2 = (next.delta>>8);
    next.delta = d1;
    backtrackb(next,T,Tb,Tm,bps);
    next.delta = d2;
    next.j=next.i-1;
    next.i=curr.i+1;
    backtrackm(next,T,Tb,Tm,bps);
  }
}

void backtrackm(struct index curr, struct index ***T,
		struct index ***Tb, struct index ***Tm, int *bps) {
  struct index next = Tm[curr.i][curr.j][curr.delta];
  int d1,d2;
  if (next.type==1) {
    backtrackb(next,T,Tb,Tm,bps);
  }
  if (next.type==2) {
    d1 = next.delta & 0xff, d2 = (next.delta>>8);
    next.delta = d1;
    backtrackb(next,T,Tb,Tm,bps);
    next.delta = d2;
    next.j = next.i-1;
    next.i = curr.i;
    backtrackm(next,T,Tb,Tm,bps);
  }
  else if (next.type==5) {
    backtrackm(next,T,Tb,Tm,bps);
  }
}

void pf(char *a) {
  int i,j,d;
  int n = strlen(a);
  int k, l;
  double **Z = (double **) xcalloc(n+1,sizeof(double *));
  double **Zb = (double **) xcalloc(n+1,sizeof(double *));
  double **Zm = (double **) xcalloc(n+1,sizeof(double *));
  double RT = 0.0019872370936902486 * (37 + 273.15) * 100; // 0.01 * (kcal K)/mol
  double energy;
  double expE;
  double tempZ;
  int IJ, JI;

  static int PN[5][5];
  int *seq = (int *)xcalloc(n+1,sizeof(int));

  const char *fname = "energy.par";
  read_parameter_file(fname);

  N = n;

  initialize_PN(PN);
  translateseq(a,seq);

  for (j=0;j<=N;j++) {
    Z[j] = (double *) xcalloc(N+1,sizeof(double));
    Zb[j] = (double *) xcalloc(N+1,sizeof(double));
    Zm[j] = (double *) xcalloc(N+1,sizeof(double));
  }

  /* Initialize Sequences of lengths shorter than or equal to THRESHOLD=3
   * will have no structure */

  for (d = 0; d <= THRESHOLD; d++) {
    for (i = 1; i <= n - d; i++) {
      j = i + d;
      Z[i][j] = 1.0;
      Zb[i][j] = 0.0;
      Zm[i][j] = 0.0;
    }
  }

  for (d = THRESHOLD+1; d <= n-1; d++) {
    for (i = 1; i <= n - d; i++) {
      j = i + d;
      Z[i][j] = 0.0;
      Zb[i][j] = 0.0;
      Zm[i][j] = 0.0;
      IJ = PN[seq[i]][seq[j]];
      JI = PN[seq[j]][seq[i]];

      /* Start by filling in Zb, this might be needed for Zm or Z */
      if (IJ) {
	tempZ = 0.0;
	/* Hairpin loop */
	energy = hairpinloop(i, j, IJ, seq, a);
	tempZ += exp(-energy/RT);
      
	for (k = i+1; k < j - THRESHOLD ; k++)
	  for (l = k + THRESHOLD + 1; l < j; l++) {
	    if (PN[seq[k]][seq[l]]) {
	      /* Interior loop, bulge or stack? */
	      energy = interiorloop(i, j, k, l, IJ, 
				    PN[seq[l]][seq[k]], seq);
	      tempZ += Zb[k][l] * exp(-energy/RT);
	      
	      /* Multi loop, where (i,j) is the closing base pair. */
	      if (k>i+THRESHOLD+2) {
		  /* If there is a multiloop with (i,j) as the closing
		   * base pair and (k,l) as the right-most base pair
		   * inside i..j, then there will be at least one
		   * 'hairpin' in the region i+1..k-1, or else it is
		   * not a multloop, but rather an interior loop. 
		   */
		energy = multiloop_closing(i,j,k,l,JI,
				   PN[seq[k]][seq[l]],seq);
		tempZ += Zb[k][l] * Zm[i+1][k-1] *
		  exp(-energy/RT);
	      }
	    }
	  }
	Zb[i][j] = tempZ;
      }

      /* Multi loop, Zm */
      tempZ = 0.0;
      for (k = i; k < j - THRESHOLD ; k++)
	for (l = k + THRESHOLD + 1; l <= j; l++) {
	  if (PN[seq[k]][seq[l]]) {
	    /* One stem */
	    energy = P->MLintern[PN[seq[k]][seq[l]]] + P->MLbase*(k-i+j-l);
	    /* Dangles */
	    if (DANGLE) {
	      if (k>1)
		/* Dangle between (k,l) and k-1 */
		energy += P->dangle5[PN[seq[k]][seq[l]]][seq[k-1]];
	      if (l<N)
		/* Dangle between (k,l) and l+1 */
		energy += P->dangle3[PN[seq[k]][seq[l]]][seq[l+1]];
	    }
	    /* terminal base pair is not GC or CG */
	    //if (PN[seq[k]][seq[l]]>2)
	    //  energy += TerminalAU;
	    tempZ += Zb[k][l] * exp(-energy/RT);
	    
	    /* More than one stem */
	    if (k>i+THRESHOLD+2) {
	      energy = P->MLintern[PN[seq[k]][seq[l]]] + P->MLbase*(j-l);
	      /* Dangles */
	      if (DANGLE) {
		/* Between (k,l) and k-1 */
		energy += P->dangle5[PN[seq[k]][seq[l]]][seq[k-1]];
		/* Between (k,l) and l+1 */
		energy += P->dangle3[PN[seq[k]][seq[l]]][seq[l+1]];
	      }
	      tempZ += Zb[k][l] * Zm[i][k-1] * exp(-energy/RT);
	    }
	  }
	}
      Zm[i][j] = tempZ;
      
      /* Z */
      /* The unpaired structure has energy 0 */
      tempZ = 1.0;
      for (k = i; k < j - THRESHOLD ; k++)
	for (l = k + THRESHOLD + 1; l <= j; l++) {
	  if (PN[seq[k]][seq[l]]) {
	    /* Compute Z */
	    energy = 0.0; 
	    /* dangle */
	    if (DANGLE) {
	      if (k>1)
		energy += P->dangle5[PN[seq[k]][seq[l]]][seq[k-1]];
	      if (l<n)
		energy += P->dangle3[PN[seq[k]][seq[l]]][seq[l+1]];
	    }
	    /* Terminal AU */
	    if (PN[seq[k]][seq[l]]>2)
	    energy += TerminalAU;
	    expE = exp(-energy/RT);
	    if (k==i) {
	      tempZ += Zb[k][l] * expE;
	    }
	    else
	      tempZ += Z[i][k-1] * Zb[k][l] * expE;
	  } /* if (PN[seq[k]][seq[l]]) */
	} /* for */
      Z[i][j] = tempZ;
    }
  } 
  printf("%.12g\n",Z[1][n]);//,-RT/100.*log(Z[1][n]));
  
  for (j=0;j<=N;j++) {
    free(Z[j]);
    free(Zb[j]);
    free(Zm[j]);
    }
  free(Z);
  free(Zb);
  free(Zm); 
}

void sub_structure(int i, int j, int *bps, int *sub_bps) {
  int k;
  sub_bps[0]=0;
  for (k = 1; k<i; k++)
    sub_bps[k] = -1;
  for (k = j+1; k<=N; k++)
    sub_bps[k] = -1;
  for (k = i; k<=j; k++) {
    if ( k<bps[k] && bps[k]<=j) {
      sub_bps[k] = bps[k];
      sub_bps[bps[k]] = k;
      sub_bps[0]++;
    }
  }
}

/* Number of base pairs in the region i to j in bps */
int numbp(int i, int j, int *bps) {
  int n=0;
  int k;
  for (k=i;k<=j;k++)
    if ( k<bps[k] && bps[k]<=j )
      n++;
  return n;
}

int bn(char A) {
/* Return the integer corresponding to the given base
   @ = 0  A = 1  C = 2  G = 3  U = 4 */
  if (A=='A')
    return 1;
  else if (A=='C')
    return 2;
  else if (A=='G')
    return 3;
  else if (A=='U')
    return 4;
  else
    return 0;
}

/* A function for printing all base pairs in a base pair list. This
 * can be used for debugging. */
void print_bps(int *bps) {
  int k;
  printf("basepairs:\n");
  for (k = 1;k<=N;k++)
    if ( k<bps[k] )
      printf("(%u,%u)\n",k,bps[k]);
  printf("\n");
}

void printbpsstring(int *bps){
  /* Print the secondary structure as a string */
  int i;
  for (i=1;i<=N;i++)
    if (bps[i]==0)
      printf(".");
    else if (i<bps[i])
      printf("(");
    else
      printf(")");
  //printf("\n");
}

void initialize_PN(int pn[5][5]) {
  int i,j;
  for (i=0;i<5;i++)
    for (j=0;j<5;j++)
      pn[i][j] = 0;
  pn[1][4]=5;
  pn[4][1]=6;
  pn[2][3]=1;
  pn[3][2]=2;
  pn[3][4]=3;
  pn[4][3]=4;
    /*
    {{0,0,0,0,0},
    {0,0,0,0,5},
    {0,0,0,1,0},
    {0,0,2,0,3},
    {0,6,0,4,0}}
    */
}

void translateseq(char *a, int *seq) {
  int i;
  seq[0] = strlen(a);
  for (i=0;i<seq[0];i++)
    seq[i+1] = bn(a[i]);
}

void initialize_NumBP(int **NumBP, int *bps, int n){
  int d,i,j;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      NumBP[i][j] = 0;
  for (d = THRESHOLD+1; d < n; d++) {
    for (i = 1; i <= n - d; i++) {
      j = i + d;
      NumBP[i][j] = numbp(i,j,bps);
    }
  }
}

double hairpinloop(int i, int j, int bp_type, int *seq, char *a) {
  double energy;
  energy = ((j-i-1) <= 30) ? P->hairpin[j-i-1] :
    P->hairpin[30]+(P->lxc*log((j-i-1)/30.));
  if ((j-i-1) >3) /* No mismatch for triloops */
    energy += P->mismatchH[bp_type]
      [seq[i+1]][seq[j-1]];
  else {
    char tl[6]={0}, *ts;
    strncpy(tl, a+i-1, 5);
    if ((ts=strstr(Triloops, tl))) 
      energy += P->Triloop_E[(ts - Triloops)/6];
    if (bp_type>2)
      energy += TerminalAU;
	}
  if ((j-i-1) == 4) { /* check for tetraloop bonus */
    char tl[7]={0}, *ts;
    strncpy(tl, a+i-1, 6);
    if ((ts=strstr(Tetraloops, tl))) 
      energy += P->TETRA_ENERGY[(ts - Tetraloops)/7];
  }
  return energy;
}

double interiorloop(int i, int j, int k, int l, int bp_type1, int bp_type2, int *seq) {
  double energy;
  int n1,n2;
  /* Interior loop, bulge or stack? */
  n1 = k-i-1;
  n2 = j-l-1;

  if ( (n1>0) && (n2>0) )
    /* Interior loop, special cases for interior loops of 
     * sizes 1+1, 1+2, 2+1, and 2+2. */
    if ( (n1==1) && (n2==1) )
      energy = P->int11[bp_type1][bp_type2]
	[seq[i+1]][seq[j-1]];
    else if ( (n1==2) && (n2==1) )
      energy = P->int21[bp_type2][bp_type1]
	[seq[j-1]][seq[i+1]][seq[i+2]];
    else if ( (n1==1) && (n2==2) )
      energy = P->int21[bp_type1][bp_type2]
	[seq[i+1]][seq[l+1]][seq[l+2]];
    else if ( (n1==2) && (n2==2) )
      energy = P->int22[bp_type1][bp_type2]
	  [seq[i+1]][seq[i+2]][seq[l+1]][seq[l+2]];
    else
      energy = ((n1+n2<=30)?(P->internal_loop[n1+n2]):
	(P->internal_loop[30]+(P->lxc*log((n1+n2)/30.))))+
	P->mismatchI[bp_type1][seq[i+1]][seq[j-1]] +
	P->mismatchI[bp_type2][seq[l+1]][seq[k-1]] +
	min2(MAX_NINIO, abs(n1-n2)*P->F_ninio[2]);
  else if ( (n1>0) || (n2>0) ) {
    /* Bulge */
    energy = ((n2+n1)<=30)?P->bulge[n1+n2]:
      (P->bulge[30]+(P->lxc*log((n1+n2)/30.)));
    if ( (n1+n2)==1 )
      /* A bulge of size one is so small that the base pairs
       * can stack */
      energy += P->stack[bp_type1][bp_type2];
    else {
      if ( bp_type1>2)
	energy += TerminalAU;
      if ( bp_type2>2)
	energy += TerminalAU;
    }
  }
  else { /* n1=n2=0 */
    /* Stacking base pair */
    energy = P->stack[bp_type1][bp_type2];
  }
  return energy;
}

double multiloop_closing(int i, int j, int k, int l, int bp_type1, int bp_type2, int *seq) {
  /* Multiloop where (i,j) is the closing base pair */
  double energy;
  energy = P->MLclosing + P->MLintern[bp_type1] + P->MLintern[bp_type2] + P->MLbase*(j-l-1);
  /* Dangles */
  if (DANGLE) {
    /* Dangles might be counted twice and also a dangle
     * between an outermost base pair might be counted also 
     * if the neighbouring base is paired. This is a 
     * simplified way to take dangles into account that is
     * also used in RNAfold -p0 
     */
    /* (i,j) to i+1 */
    energy += P->dangle3[bp_type1][seq[i+1]];
    /* (k,l) to k-1 */
    energy += P->dangle5[bp_type2][seq[k-1]];
    /* (i,j) to j-1 */  
    energy += P->dangle5[bp_type1][seq[j-1]];
    /* (k,l) to l+1 */  
    energy += P->dangle3[bp_type2][seq[l+1]];
  }
  return energy;
}

// ****************************************************************************
// FFTbor functions
// ****************************************************************************

void solveSystem(int sequenceLength, dcomplex **rootsOfUnity, double *coefficients, double scalingFactor) {
  int i;
  dcomplex sum = ZERO_C;
  
	if (FFTBOR_DEBUG) {
		printMatrix(rootsOfUnity, (char *)"START ROOTS AND SOLUTIONS", 0, sequenceLength, 0, 1);
	  std::cout << "END ROOTS AND SOLUTIONS" << std::endl << std::endl;
	}

  fftw_complex signal[sequenceLength + 1];
  fftw_complex result[sequenceLength + 1];
  
  for (i = 0; i <= sequenceLength; i++) {
    signal[i][FFTW_REAL] = (pow(10, PRECISION) * rootsOfUnity[i][1].real()) / scalingFactor;
    signal[i][FFTW_IMAG] = (pow(10, PRECISION) * rootsOfUnity[i][1].imag()) / scalingFactor;
  }
  
  fftw_plan plan = fftw_plan_dft_1d(sequenceLength + 1, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  for (i = 0; i <= sequenceLength; i++) {
    coefficients[i] = PRECISION == 0 ? result[i][FFTW_REAL] / (sequenceLength + 1) : pow(10.0, -PRECISION) * static_cast<int>(result[i][FFTW_REAL] / (sequenceLength + 1));
    sum            += coefficients[i];
    
    std::cout << i << "\t" << coefficients[i] << std::endl;
  }
  
	std::cout << "Scaling factor (Z{1, n}): " << scalingFactor << std::endl;
  std::cout << "Sum: " << sum << std::endl;
}

/* Number of base pairs in the region i to j in basePairs */
int numberOfBasePairs(int i, int j, int *basePairs) {
  int n = 0;
  int k;
  for (k = i; k <= j; ++k) {
    // If position k opens a b.p. '(' and is closed within the range [i, j], increment n
    if (k < basePairs[k] && basePairs[k] <= j) {
      n++;
    }
  }
  
  return n;
}

int** fillBasePairCounts(int *basePairs, int n) {
  int i, d, **bpCounts;
  
  bpCounts = (int **) calloc(n + 1, sizeof(int *));
  for (i = 1; i <= n; ++i) {
    bpCounts[i] = (int *) calloc(n + 1, sizeof(int));
  }
  
  for (d = MIN_PAIR_DIST + 1; d < n; ++d) {
    for (i = 1; i <= n - d; ++i) {
      bpCounts[i][i + d] = numberOfBasePairs(i, i + d, basePairs);
    }
  }
  
  return bpCounts;
}

int jPairedTo(int i, int j, int *basePairs) {
  return basePairs[i] == j ? -1 : 1;
}

int jPairedIn(int i, int j, int *basePairs) {
  return basePairs[j] >= i && basePairs[j] < j ? 1 : 0;
}

/* Checks if the i-th and j-th positions in the sequence (1-indexed!) can base pair */
int canBasePair(int i, int j, char *sequence) {	
  if (j - i <= MIN_PAIR_DIST)
    return 0;
  else if ((sequence[i] == 'A' && sequence[j] == 'U') || (sequence[i] == 'U' && sequence[j] == 'A'))
    return 1;
  else if ((sequence[i] == 'C' && sequence[j] == 'G') || (sequence[i] == 'G' && sequence[j] == 'C'))
    return 1;
  else if ((sequence[i] == 'U' && sequence[j] == 'G') || (sequence[i] == 'G' && sequence[j] == 'U'))
    return 1;
  else
    return 0;
}

void printMatrix(dcomplex **matrix, char *title, int iStart, int iStop, int jStart, int jStop) {
  int i, j;
  
  printf("%s\n", title);
      
  for (i = iStart; i <= iStop; ++i) {
    for (j = jStart; j <= jStop; ++j) {
      printf("%+.15f, %-+25.15f", matrix[i][j].real(), matrix[i][j].imag());
    }
    std::cout << std::endl;
  }
}
