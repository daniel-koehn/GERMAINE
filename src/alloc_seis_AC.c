/*
 * allocate memory for FD seismograms 
 *
 * Daniel Koehn
 * Kiel, 21/06/2016
 */

#include "fd.h"

void alloc_seis_AC(struct waveAC *waveAC, int ntr){

        /* global variables */
	extern int NX, NY, NONZERO, NXNY;

	/* local variables */
        int i;

        /* allocate memory for receiver vector */
	(*waveAC).precr = vector(1,ntr);
	(*waveAC).preci = vector(1,ntr);

}



