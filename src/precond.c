/*------------------------------------------------------------------------
 * Preconditioning of the gradient 
 * 
 * Daniel Koehn
 * Kiel, 14.12.2015
 * ----------------------------------------------------------------------*/

#include "fd.h"

void precond(float ** grad, int nsrc, float ** srcpos, int ** recpos, int ntr, int iter){

        /* global variables */
	extern int SPATFILTER;
	
        /* local variables */
	int i, j;
		
	/* Smooth gradient */
	/* --------------- */
       
	/* apply wavenumber damping */
	/*if(SPATFILTER==1){
	  wavenumber(grad);
	}*/

	/* apply smooth2 */
	if(SPATFILTER==2){
	  smooth2(grad);
	}

	/* apply median filter */
	if(SPATFILTER==3){
	  smooth_grad(grad);
	}

	/* apply Gaussian filter */
	if(SPATFILTER==4){
	  gauss_filt(grad);
	}   

	/* apply median filter at source positions */
	/*median_src(Hgrad,taper_coeff,srcpos,nsrc,recpos,ntr,iter,0);*/

}
