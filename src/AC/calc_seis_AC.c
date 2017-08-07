/*------------------------------------------------------------------------
 *  Estimate FD wavefield at receiver positions                            
 *  
 *  D. Koehn
 *  Kiel, 21/06/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void calc_seis_AC(struct waveAC *waveAC, int ** recpos, int ntr, int ishot, int nshots, int nfreq){

	/* global variables */
	extern int NPML, FSSHIFT, NF;	

	/* local variables */
	int i, index;

	for(i=1;i<=ntr;i++){
	    index = i + ntr * (ishot-1) + ntr * nshots * (nfreq-1);
	    (*waveAC).precr[index] = (*waveAC).pr[recpos[2][i]+FSSHIFT][recpos[1][i]+NPML];
            (*waveAC).preci[index] = (*waveAC).pi[recpos[2][i]+FSSHIFT][recpos[1][i]+NPML];
	}

}
