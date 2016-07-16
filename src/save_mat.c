/*
 * Store (n x m) matrix A in (n x m) matrix B
 *
 * Daniel Koehn
 * Kiel, 10/12/2015
 */

#include "fd.h"

void store_mat(float ** A, float ** B, int n, int m){

	/* local variables */
	int i, j;
	
	/* store matrix A in matrix B */
	for (i=1;i<=n;i++){
	    for (j=1;j<=m;j++){
	            
	       B[j][i] = A[j][i];
	            
	    }
	}

}



