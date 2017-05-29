/*
 * allocate memory for material parameters in elastic SH-forward problem 
 *
 * Daniel Koehn
 * Kiel, 26/05/2017
 */

#include "fd.h"

void alloc_matSH(struct matSH *matSH){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j;

        /* vs, ivs2, rho and k2 */
	(*matSH).vs = matrix(1,NY,1,NX);
	(*matSH).ivs2 = matrix(1,NY,1,NX);
	(*matSH).rho = matrix(1,NY,1,NX);
	(*matSH).ivs2 = matrix(1,NY,1,NX);
	(*matSH).k2 = matrix(1,NY,1,NX);

}
