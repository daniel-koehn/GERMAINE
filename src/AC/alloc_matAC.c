/*
 * allocate memory for material parameters in acoustic forward problem 
 *
 * Daniel Koehn
 * Kiel, 18/06/2016
 */

#include "fd.h"

void alloc_matAC(struct matAC *matAC){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j;

        /* vp, ivp2 and k2 */
	(*matAC).vp = matrix(1,NY,1,NX);
	(*matAC).ivp2 = matrix(1,NY,1,NX);
	(*matAC).k2 = matrix(1,NY,1,NX);

}
