/*
 * Calculate L1 norm of TT-TTold
 *
 * Daniel Koehn
 * Kiel, 10/12/2015
 */

#include "fd.h"

float norm1(float ** TT, float ** TTold){

	/* global variables */
	extern int NX, NY;
	
	/* local variables */
	int i, j;
        float sum, norml1;
	
	/* estimate L1 norm */
        sum = 0.0;
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){
	            
	       sum += (TT[j][i]-TTold[j][i]);
	            
	    }
	}

        norml1 = sum/(NX*NY);

    return norml1;
}



