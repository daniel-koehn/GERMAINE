/*---------------------------------------------------------------------------
 * Rotate models and gradients of l-BFGS vectors if LBFGS_pointer > NLBFGS 
 * 
 * Daniel Koehn
 * Kiel, 15.09.2016
 *
 * --------------------------------------------------------------------------*/

#include "fd.h"

void rot_LBFGS_vec(float * y_LBFGS, float * s_LBFGS, int NLBFGS, int NLBFGS_vec){
	
	int h, k;

	/* --------------------- */
	/* rotate l-BFGS vectors */
	/* --------------------- */

	h = NLBFGS_vec + 1;
	for (k=1;k<=((NLBFGS-1)*NLBFGS_vec);k++){
   	  
            y_LBFGS[k] = y_LBFGS[h];
            s_LBFGS[k] = s_LBFGS[h];  
          
	    h++;
          
       }     
	
}
