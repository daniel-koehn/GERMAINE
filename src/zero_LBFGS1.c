/*------------------------------------------------------------------------
 *   zero s- and y-LBFGS-vectors
 *  
 *  
 *   last update 14.12.2015, D. Koehn
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_LBFGS1(int NLBFGS, int NLBFGS_vec, float * y_LBFGS, float * s_LBFGS){


	register int i;

	for (i=1;i<=(NLBFGS_vec*NLBFGS);i++){
            y_LBFGS[i]=0.0;
            s_LBFGS[i]=0.0;
	}	

}
