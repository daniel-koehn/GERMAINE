/*------------------------------------------------------------------------
 *  Assemble source vector
 *
 *  D. Koehn
 *  Kiel, 20.06.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void RHS_source_AC(struct waveAC *waveAC, float ** srcpos, int ishot){

	extern int NX, NY, NXNY, NPML;
        extern float DH;
	
	/* local variables */
	int i, j, k_src, nxsrc, nysrc;
        float amp;

	/* calculate discrete source position */
	nxsrc=iround(srcpos[1][ishot]/DH);
	nysrc=iround(srcpos[2][ishot]/DH); 

	/* shift discrete source coordinates by NPML */
	nxsrc += NPML;
	nysrc += NPML;
	/*printf("nxsrc = %d \t nysrc = %d \n",nxsrc,nysrc);*/

	/* calculate source vector index k_src for shot no. ishot */
	k_src = nxsrc + (nysrc-1) * NX;

	//printf("ksrc = %d \n",k_src);

        /* spike wavelet*/
        amp = 1.0;

        /* initialize source vector */
        for (i=0;i<NXNY;i++){
	     (*waveAC).RHSr[i] = 0.0;
	     (*waveAC).RHSi[i] = 0.0;
	}

        /* Define spike source */
        (*waveAC).RHSr[k_src] = amp/(DH*DH);
        (*waveAC).RHSi[k_src] = 0.0;


}




