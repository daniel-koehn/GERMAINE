/*
 * free memory for FD seismograms 
 *
 * Daniel Koehn
 * Kiel, 21/06/2016
 */

#include "fd.h"

void free_seis_AC(struct waveAC *waveAC, int ntr){

        /* global variables */

	/* local variables */

        /* allocate memory for receiver vector */
	free_vector((*waveAC).precr, 1, ntr);
	free_vector((*waveAC).preci, 1, ntr);
}



