/*------------------------------------------------------------------------
 * Apply taper function  
 * 
 * Daniel Koehn
 * Kiel, 03.07.2016
 * ----------------------------------------------------------------------*/

#include "fd.h"

void taper_grad_hor(float ** grad){

        /* global variables */
	extern int NX, NY, GRADT2, NPML;
	extern int SWS_TAPER_GRAD_HOR;
	extern float  EXP_TAPER_GRAD_HOR, DH;
	
        /* local variables */
	int i, j;
	
	if (SWS_TAPER_GRAD_HOR){
	
	   for (j=1;j<=NY;j++){
                for (i=1;i<=NX;i++){

		    //grad[j][i] *= pow((float)(j*DH),EXP_TAPER_GRAD_HOR);

		    if(j <= GRADT2+NPML){
		       grad[j][i] = 0.0;
		    }

		}
	   }

	}

}
