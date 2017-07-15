/*------------------------------------------------------------------------
 *  Assemble source vector
 *
 *  D. Koehn
 *  Kiel, 20.06.2016
 *  last update: 15.07.2017 
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void RHS_source_AC(struct waveAC *waveAC, float ** srcpos, int ishot){

	extern int NX, NY, NXNY, NPML, FSSHIFT, READ_WAVELET;
        extern float DH;
	
	/* local variables */
	int i, j, k_src, nxsrc, nysrc;
        float *amp;

	amp = vector(1,2);

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
        amp[1] = 1.0;  /* real part */
        amp[2] = 0.0;  /* imaginary part */

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




