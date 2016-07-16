/*------------------------------------------------------------------------
 *   Output of model after end of stage
 *   
 *   Daniel Koehn
 *   last update 16.12.2015
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void model_out(float ** Vp, int iter){

	/* global variables */
	extern char INV_MODELFILE[STRING_SIZE];
	
	/* local variables */
        char modfile[STRING_SIZE];

        sprintf(modfile,"%s_vp_stage_%d.bin",INV_MODELFILE,iter);
	writemod_true(modfile,Vp,3);                                                                                                                      
                          
}

