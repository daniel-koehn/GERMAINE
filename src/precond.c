/*------------------------------------------------------------------------
 * Preconditioning of the gradient 
 * 
 * Daniel Koehn
 * Kiel, 14.12.2015
 * ----------------------------------------------------------------------*/

#include "fd.h"

void precond(float ** grad, int nsrc, float ** srcpos, int ** recpos, int ntr, int iter){

        /* global variables */
	extern int SPATFILTER, GRAD_FILTER, NX, NY;
	extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
	
        /* local variables */
	int i, j;
		
	/* Preconditioning of the gradient */
	/* ------------------------------- */
       
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

	/* apply taper functions to gradient */
	/* --------------------------------- */

	// if (SWS_TAPER_GRAD_VERT){   /*vertical gradient taper is applied*/
	//    taper_grad(grad,srcpos,nsrc,recpos,ntr,1);}

	// if (SWS_TAPER_GRAD_HOR){   /* horizontal gradient taper is applied */
	//    taper_grad(grad,srcpos,nsrc,recpos,ntr,2);
	// }

	// if (SWS_TAPER_GRAD_SOURCES){   /*cylindrical taper around sources is applied*/
	//    taper_grad(grad,srcpos,nsrc,recpos,ntr,3);}

	// if (SWS_TAPER_FILE){   /* read taper from BIN-File*/
	//    taper_grad(grad,srcpos,nsrc,recpos,ntr,4);}   

	/* apply median filter at source positions */
	/*median_src(Hgrad,taper_coeff,srcpos,nsrc,recpos,ntr,iter,0);*/

}
