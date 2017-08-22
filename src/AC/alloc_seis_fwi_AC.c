/*
 * allocate memory for field data FD seismograms 
 *
 * Daniel Koehn
 * Kiel, 07/08/2017
 */

#include "fd.h"

void alloc_seis_fwi_AC(struct waveAC *waveAC, int ntr, int nshots){

        /* global variables */
	extern int NX, NY, NONZERO, NXNY, NF;

	/* local variables */
        int i;

        /* allocate memory for receiver vector */
	(*waveAC).pobsr = vector(1,ntr*NF*nshots);
	(*waveAC).pobsi = vector(1,ntr*NF*nshots);

}



