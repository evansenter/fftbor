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
/* Vienna RNA h-files */
#include "energy_const.h"
#include "energy_par.h"
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

/* Some Vienna RNA things */
extern void  read_parameter_file(const char fname[]);

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
  int i,j,d,delta;
  int n = strlen(a);
  int k, l, e, f;
  int bpdiff;
  double ***Z = (double ***) xcalloc(n+1,sizeof(double **));
  double ***Zb = (double ***) xcalloc(n+1,sizeof(double **));
  double ***Zm = (double ***) xcalloc(n+1,sizeof(double **));
  /* Backtrack */
  struct index ***T = (struct index ***) xcalloc(n+1,sizeof(struct index **));
  struct index ***Tb = (struct index ***) xcalloc(n+1,sizeof(struct index **));
  struct index ***Tm = (struct index ***) xcalloc(n+1,sizeof(struct index **));
  double ***NB = (double ***) xcalloc(n+1,sizeof(double **));
  double RT = 0.0019872370936902486 * (temperature + 273.15) * 100; // 0.01 * (kcal K)/mol
  double energy;
  double expE = 0.0;
  double tempZ = 0.0, tempN, tempMFE;
  struct index tempT, tempt;
  double temp1, temp2;
  int IJ, JI;
  int zero_count;

  int *seq = xcalloc(n+1,sizeof(int));
  static int PN[5][5];
  static int **NumBP;

  struct index curr;
  int *bps_delta;

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

  N = n;
  for (i=0;i<=N;i++) {
    Z[i] = (double **) xcalloc(N+1,sizeof(double *));
    Zb[i] = (double **) xcalloc(N+1,sizeof(double *));
    Zm[i] = (double **) xcalloc(N+1,sizeof(double *));
    T[i] = (struct index **) xcalloc(N+1,sizeof(struct index *));
    Tb[i] = (struct index **) xcalloc(N+1,sizeof(struct index *));
    Tm[i] = (struct index **) xcalloc(N+1,sizeof(struct index *));
    if (NUMBER)
      NB[i] = (double **) xcalloc(N+1,sizeof(double *));
    for (j=0;j<=N;j++) {
      Z[i][j] = (double *) xcalloc(DELTA+1,sizeof(double));
      Zb[i][j] = (double *) xcalloc(DELTA+1,sizeof(double));
      Zm[i][j] = (double *) xcalloc(DELTA+1,sizeof(double));
      if (j>=i) {
	T[i][j] = (struct index *) xcalloc(DELTA+1,sizeof(struct index));
	Tb[i][j] = (struct index *) xcalloc(DELTA+1,sizeof(struct index));
	Tm[i][j] = (struct index *) xcalloc(DELTA+1,sizeof(struct index));
	if (NUMBER)
	  NB[i][j] = (double *) xcalloc(DELTA+1,sizeof(double));
      }
    }
  }


  /* Backtrack */
  bps_delta = (int *) xcalloc(N+1, sizeof(int));
  curr.type = 0;
  curr.i = 1;
  curr.j = N;

 /* Initialize Sequences of lengths shorter than or equal to THRESHOLD=3
  * will have no delta-neighbours, except 0-neighbours */
  for (d = 0; d <= THRESHOLD; d++) {
    for (i = 1; i <= n - d; i++) {
      j = i + d;
      delta = 0;
      if (NUMBER)
	NNB(delta,i,j) = 1.0;
      ZZ(delta,i,j) = 1.0;
      ZZb(delta,i,j) = 0.0;
      ZZm(delta,i,j) = 0.0;
#ifdef COMPUTEMFE
      MFEMFE(delta,i,j) = 0.0;
      MFEMFEb(delta,i,j) = INF;
      MFEMFEm(delta,i,j) = INF;
#endif
      T[i][j][delta].delta = 0; T[i][j][delta].i = 0; T[i][j][delta].j = 0;
      Tb[i][j][delta].delta = 0; Tb[i][j][delta].i = 0; Tb[i][j][delta].j = 0;
      Tm[i][j][delta].delta = 0; Tm[i][j][delta].i = 0; Tm[i][j][delta].j = 0;
      for (delta = 1; delta <= DELTA; delta++) {
	if (NUMBER)
	  NNB(delta,i,j) = 0.0;
	ZZ(delta,i,j) = 0.0;
	ZZb(delta,i,j) = 0.0;
	ZZm(delta,i,j) = 0.0;
#ifdef COMPUTEMFE
	MFEMFE(delta,i,j) = INF;
	MFEMFEb(delta,i,j) = INF;
	MFEMFEm(delta,i,j) = INF;
#endif
	T[i][j][delta].type = Tb[i][j][delta].type = Tm[i][j][delta].type = 0;
      }
    }
  }

  zero_count = 0; /* Count number of zeros in a row in Z[][1][n]*/
  
  for (delta = 0;delta <= DELTA; delta++) {
    if (zero_count>=STOP)
      break;
      
    for (d = THRESHOLD+1; d <= n-1; d++) {
      for (i = 1; i <= n - d; i++) {
	      j = i + d;
	      
	      if (NUMBER)
	        NNB(delta,i,j) = 0.0;
	        
	      ZZ(delta,i,j) = 0.0;
  	    ZZb(delta,i,j) = 0.0;
  	    ZZm(delta,i,j) = 0.0;

        #ifdef COMPUTEMFE
  	      MFEMFE(delta,i,j) = INF;
  	      MFEMFEb(delta,i,j) = INF;
  	      MFEMFEm(delta,i,j) = INF;
        #endif

	      T[i][j][delta].type = Tb[i][j][delta].type = Tm[i][j][delta].type = 0;
	      IJ = PN[seq[i]][seq[j]];
	      JI = PN[seq[j]][seq[i]];

	      /* Start by filling in Zb, this might be needed for Zm or Z */
	      /* Start by filling in MFEb, this might be needed for MFEm or MFE */
	      if (IJ) {
	        tempZ = 0.0;
	        tempMFE = INF;
	        tempT.type = 0;
	        
	        /* (i,j) forms a base pair */
	        if ( delta == (NumBP[i][j] + ((bps[i]==j) ? -1 : 1)) ) {
	          /* Hairpin loop */
	          tempMFE = hairpinloop(i,j,IJ,seq,a);
	          tempZ += exp(-tempMFE/RT);
	          tempT.type = 0;
	        }

	        for (k = i+1; k < j - THRESHOLD ; k++)
	          for (l = max2(k + THRESHOLD + 1,j-Linterior-1); l < j; l++) {
	            if (PN[seq[k]][seq[l]]) {
		            /* In the neighbour (k,l) is the right-most base pair in 
		             * the region from i+1 to j-1. */
		            bpdiff = NumBP[i][j] - NumBP[k][l] + ((bps[i]==j) ? -1:1);
		            
		            if (bpdiff <= delta) {
		              /* Interior loop, bulge or stack? */
		              energy = interiorloop(i, j, k, l, IJ, PN[seq[l]][seq[k]], seq);
		              tempZ += ZZb(delta - bpdiff,k,l) * exp(-energy/RT);
                  
                  #ifdef COMPUTEMFE
		                if (MFEMFEb(delta - bpdiff,k,l)<INF && mindouble(MFEMFEb(delta - bpdiff,k,l) + energy,&tempMFE)) {
		                  tempT.type = 1;
		                  tempT.delta = delta - bpdiff;
		                  tempT.i=k;
		                  tempT.j=l;
		                }
                  #endif
		            }

		            /* Multi loop, where (i,j) is the closing base pair. */
		            if (k>i+THRESHOLD+2) {
		              /* If there is a multiloop with (i,j) as the closing
		               * base pair and (k,l) as the right-most base pair
		               * inside i..j, then there will be at least one
		               * 'hairpin' in the region i+1..k-1, or else it is
		               * not a multloop, but rather an interior loop. 
		               */
		              bpdiff -= NumBP[i+1][k-1];
		              
		              if (bpdiff <= delta){
		                energy = multiloop_closing(i,j,k,l,JI, PN[seq[k]][seq[l]],seq);
		                
		                for (e=0;e<=delta-bpdiff;e++) {
		                  f = delta -bpdiff -e;
		                  tempZ += ZZb(e,k,l) * ZZm(f,i+1,k-1) * exp(-energy/RT);
                    
                      #ifdef COMPUTEMFE
		                    if (MFEMFEb(e,k,l)<INF && MFEMFEm(f,i+1,k-1)<INF && mindouble(MFEMFEb(e,k,l) + MFEMFEm(f,i+1,k-1) + energy,&tempMFE)) {
			                    tempT.type = 2;
			                    tempT.delta = (f<<8) | e;
			                    tempT.i=k;
			                    tempT.j=l;
		                    }
                      #endif
		                } /* for e */
		              } /* if (bpdiff <= delta) */
		            }
	            }
	          } /* for l */
	          
	          ZZb(delta,i,j) = tempZ;
            
            #ifdef COMPUTEMFE
	            MFEMFEb(delta,i,j) = tempMFE;
            #endif

	          Tb[i][j][delta] = tempT;
	        } /* if IJ */
	
	        /* Multi loop, Zm */
	        tempZ = 0.0;
	        tempMFE = INF;
          
          /* j form no base pair */
          if ((bps[j]>=i) & (bps[j]<j))
            bpdiff = 1;
          else
            bpdiff = 0;
            
          if (delta>=bpdiff) {
	          tempZ = ZZm(delta-bpdiff,i,j-1)*exp(-P->MLbase/RT);
	          
            #ifdef COMPUTEMFE
	            if (MFEMFEm(delta-bpdiff,i,j-1) < INF) {
	              tempMFE = MFEMFEm(delta-bpdiff,i,j-1)+P->MLbase;
	              tempT.type=5; /* Godtyckligt nummer */
	              tempT.i=i;
	              tempT.j=j-1;
	              tempT.delta=delta-bpdiff;
	            }
            #endif
	        }
	        
	        for (k = i; k < j - THRESHOLD ; k++) {
	          if (PN[seq[k]][seq[j]]) {
	            /* (k,j) is the right-most base pair in the region from i to j */
	            /* One stem */
	            bpdiff = NumBP[i][j] - NumBP[k][j];

	            if (bpdiff <=delta) {
	              energy = P->MLintern[PN[seq[k]][seq[j]]] + P->MLbase*(k-i);
	              /* Dangles */
	              if (DANGLE) {
		              if (k>1)
		                /* Dangle between (k,j) and k-1 */
		                energy += P->dangle5[PN[seq[k]][seq[j]]][seq[k-1]];
		              if (j<N)
		                /* Dangle between (k,j) and j+1 */
		                energy += P->dangle3[PN[seq[k]][seq[j]]][seq[j+1]];
	              }
	              
	              tempZ += ZZb(delta - bpdiff,k,j) * exp(-energy/RT);
                
                #ifdef COMPUTEMFE
	                if (MFEMFEb(delta - bpdiff,k,j)<INF && mindouble(MFEMFEb(delta - bpdiff,k,j) + energy, &tempMFE)) {
		                tempT.type = 1;
		                tempT.delta = delta - bpdiff;
		                tempT.i=k;
		                tempT.j=j;
	                }
                #endif
	            }
	            
	            /* More than one stem */
	            if (k>i+THRESHOLD+2) {
	              /* If k is closer to i than this, there is not room
	               * for another stem */
	              bpdiff -= NumBP[i][k-1];
	              
	              if (bpdiff <= delta) {
		              energy = P->MLintern[PN[seq[k]][seq[j]]];
		              /* Dangles */
		              
		              if (DANGLE) {
		                /* Between (k,j) and k-1 */
		                energy += P->dangle5[PN[seq[k]][seq[j]]][seq[k-1]];
		                
		                /* Between (k,j) and j+1 */
		                if (j<N)
		                  energy += P->dangle3[PN[seq[k]][seq[j]]][seq[j+1]];
		              }

		              expE = exp(-energy/RT);
		              
		              for (e=0;e<=delta-bpdiff;e++) {
		                f = delta - bpdiff -e;
		                tempZ += ZZb(e,k,j) *  ZZm(f,i,k-1) * expE;
                    
                    #ifdef COMPUTEMFE
		                  if (MFEMFEb(e,k,j)<INF && MFEMFEm(f,i,k-1)<INF && mindouble(MFEMFEb(e,k,j) + MFEMFEm(f,i,k-1) + energy, &tempMFE)) {
		                    tempT.type = 2;
		                    tempT.delta = (f<<8) | e;
		                    tempT.i=k;
		                    tempT.j=j;
		                  }
                    #endif
		              } /* for e */
	              }	      
	            }
	          }
	        } /* for k */
	        
	        ZZm(delta,i,j) = tempZ;

          #ifdef COMPUTEMFE
	          MFEMFEm(delta,i,j) = tempMFE;
          #endif
	        
	        Tm[i][j][delta] = tempT;
	
	        /* Z and NB */
	        tempZ = 0.0;
	        tempMFE = INF;
	        tempN = 0.0;

	        /* j forms no base pair */
	        if ((bps[j]>=i) & (bps[j]<j))
	          bpdiff = 1;
	        else
	          bpdiff = 0;
	          
	        if (delta>=bpdiff) {
	          tempZ = ZZ(delta-bpdiff,i,j-1);
            
            #ifdef COMPUTEMFE
	            tempMFE = MFEMFE(delta-bpdiff,i,j-1);
	            tempT.type=5; /* Godtyckligt nummer */
	            tempT.i=i;
	            tempT.j=j-1;
	            tempT.delta=delta-bpdiff;
            #endif
            
	          if (NUMBER)
	            tempN = NNB(delta-bpdiff,i,j-1);
	        }
	        
	        for (k = i; k < j - THRESHOLD ; k++) {
	          if (PN[seq[k]][seq[j]]) {
	            /* (k,j) is the right-most base pair in the region from i to j*/
	            /* Compute Z */
	            if (i<k-1)
		            bpdiff = NumBP[i][j] - NumBP[i][k-1] - NumBP[k][j];
	            else
		            bpdiff = NumBP[i][j] - NumBP[k][j];
		            
	            if (bpdiff <= delta) {
		            energy = 0.0;

		            /* dangle */
		            if (DANGLE) {
		              if (k>1)
		                energy += P->dangle5[PN[seq[k]][seq[j]]][seq[k-1]];
		              if (j<n)
		                energy += P->dangle3[PN[seq[k]][seq[j]]][seq[j+1]];
		            }

		            /* terminalAU */
		            if (PN[seq[k]][seq[j]]>2)
		              energy += TerminalAU;

		            expE = exp(-energy/RT);
		            temp1 = 0.0;
		            temp2 = INF;
		            tempt.type = 0;

		            if (k==i) {
		              temp1 += ZZb(delta - bpdiff,k,j);
                  
                  #ifdef COMPUTEMFE
		                if (mindouble(MFEMFEb(delta - bpdiff,k,j), &temp2)) {
		                  tempt.type = 1;
		                  tempt.delta = delta - bpdiff;
		                  tempt.i=k;
		                  tempt.j=j;
		                }
                  #endif
		            } else {
		              for (e=0;e<=delta - bpdiff; e++) {
		                temp1 += ZZ(e,i,k-1) * ZZb(delta - bpdiff -e,k,j);
                    
                    #ifdef COMPUTEMFE
		                  if (MFEMFE(e,i,k-1)<INF && MFEMFEb(delta - bpdiff -e,k,j)<INF && mindouble(MFEMFE(e,i,k-1) + MFEMFEb(delta - bpdiff -e,k,j),&temp2)) {
		                    tempt.type = 2;
		                    tempt.delta = (e<<8)|(delta-bpdiff-e);
		                    tempt.i=k;
		                    tempt.j=j;
		                  }
                    #endif
		              }
		            }
		            
		            tempZ += temp1 * expE;
		            
		            if (temp2<INF && mindouble(temp2 + energy, &tempMFE))
		              tempT = tempt;
	            } /* if (bpdiff <= delta */ 

	            /* Compute NB */
	            if (NUMBER) {
		            if (bps[k]!=j)
		              bpdiff = NumBP[i][j] - NumBP[i][k-1] - NumBP[k+1][j-1] + 1;
		              
		            if (bpdiff <= delta) {
		              if (k==i) {
		                tempN += NNB(delta - bpdiff,k+1,j-1);
		              } else
		                for (e=0;e<=delta - bpdiff;e++) {
		                  tempN += NNB(e,i,k-1) * NNB(delta - bpdiff -e,k+1,j-1);
		                }
		            } /* if (bpdiff <= delta */
	            }
	          } /* if (PN[seq[k]][seq[l]]) */
	        } /* for */

	        ZZ(delta,i,j) = tempZ;
          
          #ifdef COMPUTEMFE
	          MFEMFE(delta,i,j) = tempMFE;
          #endif

	        T[i][j][delta] = tempT;
	        NNB(delta,i,j) = tempN;
        } /* for i */
      } /* for d */
    
    printf("%d",delta);
    printf("\t%.12g",ZZ(delta,1,N));
    
    if (NUMBER)
      printf("\t%0.12g",NNB(delta,1,n));

    #ifdef COMPUTEMFE
      if (MFEMFE(delta,1,n)<INF) {
        /* Backtrack */
        for (i=0;i<=N;i++) 
          bps_delta[i]=0;

        curr.delta = delta;
        backtrack(curr,T,Tb,Tm,bps_delta);
        printf("\t%.12g\t",MFEMFE(delta,1,n)/100);
        printbpsstring(bps_delta);
        printf("\t%.12g\n",exp(-MFEMFE(delta,1,n)/RT)/ZZ(delta,1,n));
        zero_count = 0;
      } else {
        printf("\t-\t-\n");
        zero_count++;
      }
    #else
      if (ZZ(delta,1,n)>0)
        zero_count = 0;
      else
        zero_count++;
      printf("\n");
    #endif
    
    fflush(stdout);
  }  /* for delta */

  free(bps_delta);

  for (i=1;i<=N;i++) {
    for (j=i;j<=N;j++) {
      free(Z[i][j]);
      free(Zb[i][j]);
      free(Zm[i][j]);
      free(T[i][j]);
      free(Tb[i][j]);
      free(Tm[i][j]);
      if (NUMBER)
	      free(NB[i][j]);
    }
    
    free(Z[i]);
    free(Zb[i]);
    free(Zm[i]);
    free(T[i]);
    free(Tb[i]);
    free(Tm[i]);
    if (NUMBER)
      free(NB[i]);
  }
  
  free(Z);
  free(Zb);
  free(Zm);
  free(T);
  free(Tb);
  free(Tm);
  
  if (NUMBER)
    free(NB);
  
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
  int *seq = xcalloc(n+1,sizeof(int));

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


