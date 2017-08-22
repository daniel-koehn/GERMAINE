/*
 * allocate memory for TE-mode FWI problem 
 *
 * Daniel Koehn
 * Kiel, 20/08/2017
 */

#include "fd.h"

void alloc_fwiTE(struct fwiTE *fwiTE, int ntr){

        /* global variables */
	extern int NX, NY;

	/* local variables */

      	/* memory for gradient */
   	(*fwiTE).grad_sigma = matrix(1,NY,1,NX);
   	(*fwiTE).gradm_sigma = matrix(1,NY,1,NX);

   	(*fwiTE).grad_epsilon = matrix(1,NY,1,NX);
   	(*fwiTE).gradm_epsilon = matrix(1,NY,1,NX);

   	(*fwiTE).hess_sigma = matrix(1,NY,1,NX);
   	(*fwiTE).hess_epsilon = matrix(1,NY,1,NX);

   	(*fwiTE).Hgrad_sigma = matrix(1,NY,1,NX);
   	(*fwiTE).Hgrad_epsilon = matrix(1,NY,1,NX);

   	(*fwiTE).sigma_old =  matrix(1,NY,1,NX);
   	(*fwiTE).epsilon_old =  matrix(1,NY,1,NX);

   	(*fwiTE).forwardr =  matrix(1,NY,1,NX);
   	(*fwiTE).forwardi =  matrix(1,NY,1,NX);

}



