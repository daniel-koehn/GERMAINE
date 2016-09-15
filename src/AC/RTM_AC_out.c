/*------------------------------------------------------------------------
 *   Output of RTM P-wave image
 *   
 *   Daniel Koehn
 *   last update 19.08.2016
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void RTM_AC_out(float ** Vp){

	/* global variables */
	extern char INV_MODELFILE[STRING_SIZE];
	
	/* local variables */
        char modfile[STRING_SIZE];

        sprintf(modfile,"%s_p_image.bin",INV_MODELFILE);
	writemod_true(modfile,Vp,3);                                                                                                                      
                          
}

