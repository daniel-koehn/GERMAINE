/*------------------------------------------------------------------------
 *  Calculate initial model values vp, ivp2, k2
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
		

			 (*matAC).ivp2[j][i] = 1.0/((*matAC).vp[j][i]*(*matAC).vp[j][i]);	
			 (*matAC).k2[j][i] = (*waveAC).omega2/((*matAC).vp[j][i]*(*matAC).vp[j][i]);	

			 //printf("j=%d \t i = %d \t vp = %e \n",j,i,(*matAC).ivp2[j][i]);

		}
	}

}




