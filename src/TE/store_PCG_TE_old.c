/*------------------------------------------------------------------------
 * Module for storage of Preconditioned Conjugate Gradient Method (PCG)
 * for the TE-mode case
 * 
 * Daniel Koehn
 * Kiel, 21.08.2017
 * ----------------------------------------------------------------------*/

#include "fd.h"

void store_PCG_TE_old(float * PCG_old, struct fwiTE *fwiTE){

	extern int NX, NY, IDX, IDY;
	int i, j, h;
	
        h=1;
	/* store sigma gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 PCG_old[h] = (*fwiTE).gradm_sigma[j][i];

                 h++;
	   }
	}

	/* store epsilon gradient */
	for (i=1;i<=NX;i=i+IDX){
	   for (j=1;j<=NY;j=j+IDY){

		 PCG_old[h] = (*fwiTE).gradm_epsilon[j][i];

                 h++;
	   }
	}


}
