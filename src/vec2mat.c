/*
 * Convert solution vector xr/xi to matrix pr/pi 
 *
 * Daniel Koehn
 * Kiel, 21/06/2016
 */

#include "fd.h"

void vec2mat(float **pr, float **pi, double *xr, double *xi){

        /* global variables */
	extern int NX, NY, NXNY;

	/* local variables */
        int i, j, h;

        /* convert vector to matrix */
	h = 1;
	for (j=1;j<=NY;j++){		
	   for (i=1;i<=NX;i++){

		pr[j][i] = (float) xr[h];
		pi[j][i] = (float) xi[h];
		h++;
           }
        }

}



