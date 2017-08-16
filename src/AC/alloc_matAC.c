/*
 * allocate memory for material parameters in acoustic forward problem 
 *
 * Daniel Koehn
 * Kiel, 14/08/2017
 */

#include "fd.h"

void alloc_matAC(struct matAC *matAC){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j;

        /* vp, rho and ivp2 */
	(*matAC).vp = matrix(1,NY,1,NX);
	(*matAC).rho = matrix(1,NY,1,NX);
	(*matAC).ivp2 = matrix(1,NY,1,NX);
	(*matAC).b = matrix(0,NY+1,0,NX+1);

}
