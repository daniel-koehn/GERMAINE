/*
 * Copy gradient values into boundary frame  
 *
 * Daniel Koehn
 * Kiel, 11/01/2016
 */

#include "fd.h"

void cp_grad_frame(float ** A){

        /* global variables */
	extern int NX, NY, NX0, NY0, NPML;

	/* local variables */
	int i, j, jj, PML_grad;

        /* PML_grad = 1 - copy gradient values from computation domain into PML */
	/* PML_grad = 2 - set gradient values in PML to zero */
	
	PML_grad = 2;

	/* extend model */
        /* ------------ */
	
	if(PML_grad==1){
	
	    /* left boundary */
	    for (i=1;i<=NPML;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] = A[j][NPML+1];				
		    }
	    }
	
	    /* right boundary */
	    for (i = NPML + NX0 + 1;i <= NX;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] = A[j][NPML + NX0];				
		    }
	    }

	    /* top boundary */
	    for (i=1;i<=NX;i++){
		    for (j=1;j<=NPML;j++){
			    A[j][i] = A[NPML+1][i];				
		    }
	     }

	     /* bottom boundary */
	     for (i=1;i<=NX;i++){
		    for (j = NPML + NY0 + 1;j <= NY;j++){
			    A[j][i] = A[NPML+NY0][i];				
		    }
	     }
	
	}
	
	/* set gradient in PML to zero */
	/* --------------------------- */
	
	if(PML_grad==2){
	
	    /* left boundary */
	    for (i=1;i<=NPML;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] = 0.0;				
		    }
	    }
	
	    /* right boundary */
	    for (i = NPML + NX0 + 1;i <= NX;i++){
		   for (j=1;j<=NY;j++){
			    A[j][i] = 0.0;				
		   }
	    }

	    /* top boundary */
	    for (i=1;i<=NX;i++){
		   for (j=1;j<=NPML;j++){
			   A[j][i] = 0.0;				
		   }
	    }

	    /* bottom boundary */
	    for (i=1;i<=NX;i++){
		    for (j = NPML + NY0 + 1;j <= NY;j++){
			   A[j][i] = 0.0;				
		    }
	    }
	
	}
	
			

}
