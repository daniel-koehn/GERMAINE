/*
 * allocate memory for FD seismograms 
 *
 * Daniel Koehn
 * Kiel, 21/06/2016
 */

#include "fd.h"

void alloc_seis_AC(struct waveAC *waveAC, int ntr, int nshots){

        /* global variables */
	extern int NX, NY, NONZERO, NXNY, NF;

	/* local variables */
        int i;

        /* allocate memory for receiver vector */
	(*waveAC).precr = vector(1,ntr*NF*nshots);
	(*waveAC).preci = vector(1,ntr*NF*nshots);

}



