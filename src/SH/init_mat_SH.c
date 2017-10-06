/*------------------------------------------------------------------------
 *  Calculate shear modulus
 *
 *  D. Koehn
 *  Kiel, 05.09.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void init_mat_SH(struct waveAC *waveAC, struct matSH *matSH){

	extern int NX, NY;
	
	/* local variables */
	int i, j;

	/* loop over global grid */
	for (j=1;j<=NY;j++){
	       for (i=1;i<=NX;i++){
		
			 /* calculate shear modulus */
			 (*matSH).mu[j][i] = (*matSH).rho[j][i] * (*matSH).vs[j][i] * (*matSH).vs[j][i];	

		}
	}

	/* extend buoyancy model to simplify later arithmetic averaging */
	for (j=1;j<=NY;j++){(*matSH).mu[j][0] = (*matSH).mu[j][1];}      /* left */
	for (j=1;j<=NY;j++){(*matSH).mu[j][NX+1] = (*matSH).mu[j][NX];}  /* right */

	for (i=0;i<=NX+1;i++){(*matSH).mu[0][i] = (*matSH).mu[1][i];}      /* top */
	for (i=0;i<=NX+1;i++){(*matSH).mu[NY+1][i] = (*matSH).mu[NY][i];}  /* bottom */
 
}




