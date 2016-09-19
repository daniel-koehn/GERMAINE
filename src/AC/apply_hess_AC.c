/*
 * Apply diagonal elements of Hessian
 *
 * Daniel Koehn
 * Kiel, 02.07.2016
 *
 */

#include "fd.h"

void apply_hess_AC(float ** grad, float ** hess){

        /* global variables */
	extern int NX, NY;
	extern float MAX_HESS, EPS_HESS;

	/* local variables */
	int i, j;

	/* Estimate maximum Hessian value */
        MAX_HESS = 0.0;
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){   

	       if(fabs(hess[j][i]) > MAX_HESS){MAX_HESS = fabs(hess[j][i]);}

	    }
	}
	
	/* Apply stabilized Hessian approximation to gradient */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

               grad[j][i] *= 1.0 / (hess[j][i] + EPS_HESS * MAX_HESS);    

	    }
	}
	
}



