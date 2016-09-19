/*------------------------------------------------------------------------
 *  Assemble source vector for Hessian computation
 *
 *  D. Koehn
 *  Kiel, 19.09.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void RHS_source_AC_hess(struct waveAC *waveAC, struct fwiAC *fwiAC, int ** recpos, int trace){

	extern int NX, NY, NXNY, NPML;
        extern float DH;
	
	/* local variables */
	int i, j, k_rec, nxrec, nyrec;
        float amp;

        /* initialize source vector */
        for (i=0;i<NXNY;i++){
	     (*waveAC).RHSr[i] = 0.0;
	     (*waveAC).RHSi[i] = 0.0;
	}

        /* spike wavelet*/
        amp = 1.0;

        /* calculate receiver vector index k_rec for receiver no. i */
        k_rec = (recpos[1][trace] + NPML) + (recpos[2][trace] + NPML - 1) * NX;

        /* Define spike source */
        (*waveAC).RHSr[k_rec] = amp/(DH*DH);
        (*waveAC).RHSi[k_rec] = 0.0;

}




