/*------------------------------------------------------------------------
 *  Assemble source vector for Hessian computation
 *
 *  D. Koehn
 *  Kiel, 01.09.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void RHS_source_TE_hess(struct waveAC *waveAC, int ** recpos, int trace){

	extern int NX, NY, NXNY, NPML, FSSHIFT;
        extern float DH, S;
	
	/* local variables */
	int i, j, k_rec, nxrec, nyrec;
        float ampr, ampi, mu0;
	complex float omega;

	omega = 2.0 * M_PI * ((*waveAC).freq + S * I);

        /* initialize source vector */
        for (i=0;i<NXNY;i++){
	     (*waveAC).RHSr[i] = 0.0;
	     (*waveAC).RHSi[i] = 0.0;
	}

        /* spike wavelet*/
        ampr = -creal( (omega * mu0 * I) * (1.0 + 0.0 * I) );  /* real part */
        ampi = -cimag( (omega * mu0 * I) * (1.0 + 0.0 * I) );  /* imaginary part */

        /* calculate receiver vector index k_rec for receiver no. i */
        k_rec = (recpos[1][trace] + NPML) + (recpos[2][trace] + FSSHIFT - 1) * NX;

        /* Define spike source */
        (*waveAC).RHSr[k_rec] = 1.0 / (DH*DH);
        (*waveAC).RHSi[k_rec] = 0.0 / (DH*DH);

}




