/*
 * allocate memory for field data FD seismograms 
 *
 * Daniel Koehn
 * Kiel, 07/08/2017
 */

#include "fd.h"

void alloc_seis_fwi_AC(struct fwiAC *fwiAC, int ntr, int nshots){

        /* global variables */
	extern int NX, NY, NONZERO, NXNY, NF;

	/* local variables */
        int i;

        /* allocate memory for receiver vector */
	(*fwiAC).pobsr = vector(1,ntr*NF*nshots);
	(*fwiAC).pobsi = vector(1,ntr*NF*nshots);

}



