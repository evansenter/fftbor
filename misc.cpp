/*************************************************
misc.c
P. Clote
Miscellaneous routines for reading in RNA either as
string from command line or FASTA file, etc.

NOTES:
1) I chose a[] to be an array of short int, rather than
   char, since a[0] is the length of the original RNA
   sequence, stored in a[1],...,a[ a[0] ]
   Thus this requires a special irintRNA() function.

*************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>     //for INT_MAX
#include "misc.h"
#include <ctype.h>   	// character handling, eg toupper()
#include <stdbool.h>

void *xcalloc(size_t n, size_t s) {
  void *out = calloc(n,s);
  if (!out) {
    perror("calloc error");
    exit(1);
  }
  return out;
}

int *getBasePairList(char *secStr) {
  /* Returns list L of ordered pairs (i,j) where i<j and
   * positions i,j occupied by balancing parentheses
   * For linear time efficiency, use stack
   * Assume that secStr is string consisting of '(',')' and '.'
   * Values -2,-1 returned mean NOT well balanced
   * -2 means too many ( with respect to )
   * -1 means too many ) with respect to (
   * If 1,-1 not returned, then return (possibly empty) list */
  
  int len = strlen(secStr);
  int *S = (int *) xcalloc(len/2,sizeof(int));  //empty stack
  int *L = (int *) xcalloc(2*len*(len-1)/2+1, sizeof(int)); /* initially empty
							     * list of base 
							     * pairs */
  int j, k = 0;
  char ch;

  /* First position holds the number of base pairs */
  L[0] = 0;
  for (j=1;j<=len;j++)
    L[j] = -1;

  for (j=1;j<=len;j++) {
    ch = secStr[j-1];
    if (ch == '(')
      S[k++] = j;
    else if (ch == ')') {
      if (k==0) {
	/* There is something wrong with the structure. */
	L[0] = -1; 
	return L;
      }
      else {
        L[S[--k]] = j;
	L[j] = S[k];
	L[0]++;
      }
    }
  }

  if (k != 0) {
    /* There is something wrong with the structure. */
    L[0] = -2;
  }
  
  free(S);

  return L;
}

int mindouble(double m, double *n) {
  if (m>=(*n))
    return 0;
  (*n) = m;
  return 1;
}



int min2(int m, int n) {
  if (m>n)
    return n;
  return m;
}

int max2(int m, int n) {
  if (m<n)
    return n;
  return m;
}
