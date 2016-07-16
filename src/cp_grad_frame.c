/*
 * Copy gradient values into boundary frame  
 *
 * Daniel Koehn
 * Kiel, 11/01/2016
 */

#include "fd.h"

void cp_grad_frame(float ** A){

        /* global variables */
	extern int NX, NY, NX0, NY0, NPML;

	/* local variables */
	int i, j, jj;

	/* extend model */

	/* left boundary */
	for (i=1;i<=NPML;i++){
		for (j=1;j<=NY;j++){
			A[j][i] = A[j][NPML+1];				
		}
	}
	
	/* right boundary */
	for (i = NPML + NX0 + 1;i <= NX;i++){
		for (j=1;j<=NY;j++){
			A[j][i] = A[j][NPML + NX0];				
		}
	}

	/* top boundary */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NPML;j++){
			A[j][i] = A[NPML+1][i];				
		}
	}

	/* bottom boundary */
	for (i=1;i<=NX;i++){
		for (j = NPML + NY0 + 1;j <= NY;j++){
			A[j][i] = A[NPML+NY0][i];				
		}
	}		

}
