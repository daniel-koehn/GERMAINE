/*
 * Copy sparse matrix vectors
 *
 * Daniel Koehn
 * Kiel, 11/01/2016
 */

#include "fd.h"

void cp_vec(struct waveAC *waveAC, int * Ti, int * Tj, double * Tx, double * Tz){

        /* global variables */
	extern int NONZERO;

	/* local variables */
	int i;

	/* Copy top boundary frame */
	for (i=0;i<NONZERO;i++){

		Ti[i] = (*waveAC).irow[i];	    
		Tj[i] = (*waveAC).icol[i];

		Tx[i] = (*waveAC).Ar[i];	    
		Tz[i] = (*waveAC).Ai[i];
	         
	}
	
}
