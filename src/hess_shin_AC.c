/*
 * Calculate Hessian approximation according to Shin et al. (2001) 
 *
 * Daniel Koehn
 * Kiel, 02.07.2016
 *
 */

#include "fd.h"

void hess_shin_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct matAC *matAC, float ** hess){

        /* global variables */
	extern int NX, NY, SWS_TAPER_CIRCULAR_PER_SHOT;
	extern float MAX_HESS;

	/* local variables */
	int i, j;
	float ivp6;
	
	/* assemble gradient for one shot and estimate maximum Hessian value */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

	       ivp6 = 1.0 / pow((*matAC).vp[j][i],6.0);
               hess[j][i] += 4.0 * pow((*waveAC).omega2,2.0) * creal( pow(fabs((*fwiAC).forwardr[j][i] + (*fwiAC).forwardi[j][i] * I), 4.0) * ivp6);

	    }
	}

}



