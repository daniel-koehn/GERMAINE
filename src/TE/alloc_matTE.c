/*
 * allocate memory for material parameters in TE-mode forward problem 
 *
 * Daniel Koehn
 * Kiel, 19/08/2017
 */

#include "fd.h"

void alloc_matTE(struct matTE *matTE){

        /* global variables */
	extern int NX, NY;

	/* local variables */
	int i, j;

        /* Sigma and Epsilon */
	(*matTE).sigma = matrix(1,NY,1,NX);
	(*matTE).epsilon = matrix(1,NY,1,NX);

}
