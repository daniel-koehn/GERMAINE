/*------------------------------------------------------------------------
 *  Assemble source vector for TE case
 *
 *  D. Koehn
 *  Kiel, 20.06.2016
 *  last update: 01.09.2017 
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void RHS_source_TE(struct waveAC *waveAC, float ** srcpos, int ishot){

	extern int NX, NY, NXNY, NPML, FSSHIFT, READ_WAVELET;
        extern float DH, S;
	
	/* local variables */
	int i, j, k_src, nxsrc, nysrc;
        float *amp, mu0;
	complex float omega;

	amp = vector(1,2);
	omega = 2.0 * M_PI * ((*waveAC).freq - S * I);
	mu0 = 4.0 * M_PI * 1e-7;

	/* calculate discrete source position */
	nxsrc=iround(srcpos[1][ishot]/DH);
	nysrc=iround(srcpos[2][ishot]/DH); 

	/* shift discrete source coordinates by NPML */
	nxsrc += NPML;
	nysrc += FSSHIFT;
	/*printf("nxsrc = %d \t nysrc = %d \n",nxsrc,nysrc);*/

	/* calculate source vector index k_src for shot no. ishot */
	k_src = nxsrc + (nysrc-1) * NX;

	//printf("ksrc = %d \n",k_src);

        /* spike wavelet*/
        amp[1] = -creal( (omega * mu0 * I) * (1.0 + 0.0 * I) );  /* real part */
        amp[2] = -cimag( (omega * mu0 * I) * (1.0 + 0.0 * I) );  /* imaginary part */

	/* external source wavelet */
	if(READ_WAVELET==1){
	    read_stf_dft(waveAC, amp);
	}

        /* initialize source vector */
        for (i=0;i<NXNY;i++){
	     (*waveAC).RHSr[i] = 0.0;
	     (*waveAC).RHSi[i] = 0.0;
	}

        /* Define spike source */
        (*waveAC).RHSr[k_src] = amp[1] / (DH*DH);
        (*waveAC).RHSi[k_src] = amp[2] / (DH*DH);

        /* free memory */
	free_vector(amp,1,2);

}




