/*------------------------------------------------------------------------
 *  Calculate number of non-zero elements in impedance matrix
 *
 *  D. Koehn
 *  Kiel, 18.06.2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void calc_nonzero(){

	extern int NX, NY, NONZERO;
	
	/* local variables */
	int i, j;

	/* estimate number of non-zero elements in impedance matrix */
        NONZERO = 0;
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		    /* NW gridpoint */
		    if((i > 1) && (j > 1)){NONZERO++;}
            
		    /* N gridpoint */
		    if(j > 1){NONZERO++;}
            
		    /* NE gridpoint */
		    if((i < NX) && (j > 1)){NONZERO++;}
            
		    /* W gridpoint */
		    if(i > 1){NONZERO++;}
            
		    /* central gridpoint */
		    NONZERO++;
            
		    /* E gridpoint */
		    if(i < NX){NONZERO++;}
            
		    /* SW gridpoint */
		    if( (j < NY) && (i > 1) ){NONZERO++;}
            
		    /* S gridpoint */
		    if(j < NY){NONZERO++;}
            
		    /* SE gridpoint */
		    if((i < NX) && (j < NY)) {NONZERO++;}
                       
		}
	}

}




