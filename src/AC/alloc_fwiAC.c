/*
 * allocate memory for acoustic FWI problem 
 *
 * Daniel Koehn
 * Kiel, 23/06/2016
 */

#include "fd.h"

void alloc_fwiAC(struct fwiAC *fwiAC, int ntr){

        /* global variables */
	extern int NX, NY;

	/* local variables */

      	/* memory for gradient */
   	(*fwiAC).lam =  matrix(1,NY,1,NX);
   	(*fwiAC).grad = matrix(1,NY,1,NX);
   	(*fwiAC).gradm = matrix(1,NY,1,NX);
   	(*fwiAC).hess = matrix(1,NY,1,NX);
   	(*fwiAC).Hgrad = matrix(1,NY,1,NX);
   	(*fwiAC).vp_old =  matrix(1,NY,1,NX);
   	(*fwiAC).forwardr =  matrix(1,NY,1,NX);
   	(*fwiAC).forwardi =  matrix(1,NY,1,NX);

}



