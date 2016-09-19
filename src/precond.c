/*
 * Preconditioning of gradient
 * 
 * D. Koehn
 * Kiel, 19.09.2016
 */

#include "fd.h"

void precond(float ** grad){

    /* global variables */
    extern int SPATFILTER, SWS_TAPER_GRAD_HOR, SWS_TAPER_FILE;

    /* local variables */

    /* smooth gradient with Gaussian*/
    if(SPATFILTER==1){gauss_filt(grad);}

    /* apply different tapers to gradient */
    if(SWS_TAPER_GRAD_HOR){taper_grad_hor(grad);}
    if(SWS_TAPER_FILE){taper_grad_file(grad);}
		
}

