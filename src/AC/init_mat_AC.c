/*------------------------------------------------------------------------
 *  Calculate initial model values ivp2, b
 *
 *  D. Koehn
 *  Kiel, 18.06.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void init_mat_AC(struct waveAC *waveAC, struct matAC *matAC){

	extern int NX, NY;
	
	/* local variables */
	int i, j;

	/* loop over global grid */
	for (j=1;j<=NY;j++){
	       for (i=1;i<=NX;i++){
		
			 /* calculate inverse bulk modulus */
			 (*matAC).ivp2[j][i] = 1.0 / ( (*matAC).rho[j][i] * (*matAC).vp[j][i] * (*matAC).vp[j][i] );	

			 /* calculate buoyancy */
			 (*matAC).b[j][i] = 1.0 / (*matAC).rho[j][i];

			 // printf("j=%d \t i = %d \t b = %e \n",j,i,(*matAC).b[j][i]);

		}
	}

	/* extend buoyancy model to simplify later arithmetic averaging */
	for (j=1;j<=NY;j++){(*matAC).b[j][0] = (*matAC).b[j][1];}      /* left */
	for (j=1;j<=NY;j++){(*matAC).b[j][NX+1] = (*matAC).b[j][NX];}  /* right */

	for (i=1;i<=NX;i++){(*matAC).b[0][i] = (*matAC).b[1][i];}      /* top */
	for (i=1;i<=NX;i++){(*matAC).b[NY+1][i] = (*matAC).b[NY][i];}  /* bottom */
 
}




