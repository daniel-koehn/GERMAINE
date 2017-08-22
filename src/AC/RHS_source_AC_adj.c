/*------------------------------------------------------------------------
 *  Assemble source vector for adjoint wavefield computation
 *
 *  D. Koehn
 *  Kiel, 29.06.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void RHS_source_AC_adj(struct waveAC *waveAC, int ** recpos, int ntr){

	extern int NX, NY, NXNY, NPML, FSSHIFT;
        extern float DH;
	
	/* local variables */
	int i, j, k_rec, nxrec, nyrec;
        float amp;

        /* initialize source vector */
        for (i=0;i<NXNY;i++){
	     (*waveAC).RHSr[i] = 0.0;
	     (*waveAC).RHSi[i] = 0.0;
	}

	for(i=1;i<=ntr;i++){

	    /* calculate receiver vector index k_rec for receiver no. i */
	    k_rec = (recpos[1][i] + NPML) + (recpos[2][i] + FSSHIFT - 1) * NX;

	    // printf("k_rec = %d \n",k_rec);

            /* Define adjoint source */
	    (*waveAC).RHSr[k_rec] = (*waveAC).presr[i]/(DH*DH);
            (*waveAC).RHSi[k_rec] = -(*waveAC).presi[i]/(DH*DH);

	}

	//printf("k_rec = %d \n",k_rec);

}




