/*
 * allocate memory for acoustic forward problem 
 *
 * Daniel Koehn
 * Kiel, 18/06/2016
 */

#include "fd.h"

void alloc_waveAC(struct waveAC *waveAC, struct PML_AC *PML_AC){

        /* global variables */
	extern int NX, NY, NONZERO, NXNY;

	/* local variables */

        /* impedance matrix and wavefield variables */
	(*waveAC).irow = ivector(0,NONZERO-1);
	(*waveAC).icol = ivector(0,NONZERO-1);

	(*waveAC).Ar = dvector(0,NONZERO-1);
	(*waveAC).Ai = dvector(0,NONZERO-1);

	(*waveAC).pr = matrix(1,NY,1,NX);
	(*waveAC).pi = matrix(1,NY,1,NX);

	(*waveAC).RHSr = dvector(0,NXNY-1);
	(*waveAC).RHSi = dvector(0,NXNY-1);

	(*PML_AC).Ar = matrix(1,NY,1,NX);
	(*PML_AC).Ai = matrix(1,NY,1,NX);

	(*PML_AC).Br = matrix(1,NY,1,NX);
	(*PML_AC).Bi = matrix(1,NY,1,NX);

	(*PML_AC).Cr = matrix(1,NY,1,NX);
	(*PML_AC).Ci = matrix(1,NY,1,NX);

	(*PML_AC).Axr = matrix(1,NY,1,NX);
	(*PML_AC).Axi = matrix(1,NY,1,NX);

	(*PML_AC).Byr = matrix(1,NY,1,NX);
	(*PML_AC).Byi = matrix(1,NY,1,NX);

	(*PML_AC).d_x = vector(1,NX);
	(*PML_AC).d_y = vector(1,NY);

	(*PML_AC).b_x = vector(1,NX);
	(*PML_AC).b_y = vector(1,NY);

	(*PML_AC).a_xr = vector(1,NX);
	(*PML_AC).a_yr = vector(1,NY);

	(*PML_AC).a_xi = vector(1,NX);
	(*PML_AC).a_yi = vector(1,NY);


}



