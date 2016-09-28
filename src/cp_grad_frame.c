/*
 * Copy gradient values into boundary frame  
 *
 * Daniel Koehn
 * Kiel, 11/01/2016
 */

#include "fd.h"

void cp_grad_frame(float ** A){

        /* global variables */
	extern int NX, NY, NX0, NY0, NPML, MYID, FREE_SURF, FSSHIFT;

	/* local variables */
	int i, j, jj, ii, PML_grad, yb, ye, xb, xe;
	float a, amp, DAMPING;
	char modfile[STRING_SIZE];

        /* PML_grad = 1 - copy gradient values from computation domain into PML */
	/* PML_grad = 2 - set gradient values in PML to zero */
	/* PML_grad = 3 - damp gradient values in PML with taper function */
	
	PML_grad = 3;

	/* extend model */
        /* ------------ */
	
	if(PML_grad==1){
	
	    /* left boundary */
	    for (i=1;i<=NPML;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] = A[j][NPML+1];
			    if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}				
		    }
	    }
	
	    /* right boundary */
	    for (i = NPML + NX0 + 1;i <= NX;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] = A[j][NPML + NX0];
			    if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}				
		    }
	    }

	    /* top boundary */
	    if(FREE_SURF==0){
	        for (i=1;i<=NX;i++){
		    for (j=1;j<=NPML;j++){
	                A[j][i] = A[NPML+1][i];				
			if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}
		    }
	        }
            }

	    /* bottom boundary */
	    for (i=1;i<=NX;i++){
		   for (j = FSSHIFT + NY0 + 1;j <= NY;j++){
			   A[j][i] = A[FSSHIFT+NY0][i];				
			   if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}
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
	    if(FREE_SURF==0){
	        for (i=1;i<=NX;i++){
		       for (j=1;j<=NPML;j++){
		           A[j][i] = 0.0;				
		       }
	        }
	    }

	    /* bottom boundary */
	    for (i=1;i<=NX;i++){
		    for (j = FSSHIFT + NY0 + 1;j <= NY;j++){
			   A[j][i] = 0.0;				
		    }
	    }
	
	}

	/* damp gradient in PML with taper */
	/* ------------------------------- */
	
	if(PML_grad==3){

	    float ** absorb_coeff = NULL, *coeff;

	    absorb_coeff = matrix (1,NY,1,NX);
            coeff=vector(1,NPML);

	    /* Define damping profile within PML */
	    /* --------------------------------- */ 	
	    DAMPING = 60.0;	
	    amp=1.0-DAMPING/100.0;   

    	    coeff=vector(1,NPML); 	
	    a=sqrt(-log(amp)/((NPML-1)*(NPML-1))); 		

	    for (i=1;i<=NPML;i++) 		

	        coeff[i]=exp(-(a*a*(NPML-i)*(NPML-i))); 		

		if (MYID==0){ 		
		    /*printf(" Table of coefficients \n # \t coeff \n"); 		
		    printf(" NPML=%d \t a=%f amp=%f \n", NPML,a,amp); 		

		    for (i=1;i<=NPML;i++) fprintf(FP," %d \t %5.3f \n", i, coeff[i]); */
		}			

	    /* initialize array of coefficients with one */ 	
	    for (j=1;j<=NY;j++) 	
	        for (i=1;i<=NX;i++) 
		    absorb_coeff[j][i]=1.0; 	


	    /* compute coefficients for left and right grid boundaries (x-direction) */
	    yb=1; ye=NY; 
	    for (i=1;i<=NPML;i++){
	        yb=i;
	        ye=NY-i+1;
		for (j=yb;j<=ye;j++) absorb_coeff[j][i]=coeff[i];
	    }
				    
	    yb=1; ye=NY;
	    for (i=1;i<=NPML;i++){
	        ii=NX-i+1;
		yb=i;
	        ye=NY-i+1;
		for (j=yb;j<=ye;j++) absorb_coeff[j][ii]=coeff[i];
	    }
	    	
	    /* compute coefficients for top and bottom grid boundaries (y-direction) */
	    if(FREE_SURF==0){
	        xb=1; xe=NX;
	        for (j=1;j<=NPML;j++){
	            xb=j;
		    xe=NX-j+1;
		    for (i=xb;i<=xe;i++) absorb_coeff[j][i]=coeff[j];
	        }
	    }
	    	    
	    xb=1; xe=NX;
	    for (j=1;j<=NPML;j++){
	        jj=NY-j+1;
		xb=j;
		xe=NX-j+1;
		for (i=xb;i<=xe;i++) absorb_coeff[jj][i]=coeff[j];
	    }

	    /* Apply gradient damping */
	    /* ---------------------- */

	    /* left boundary */
	    for (i=1;i<=NPML;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] *= absorb_coeff[j][i];
			    if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}				
		    }
	    }
	
	    /* right boundary */
	    for (i = NPML + NX0 + 1;i <= NX;i++){
		    for (j=1;j<=NY;j++){
			    A[j][i] *= absorb_coeff[j][i];
			    if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}				
		    }
	    }

	    /* top boundary */
	    if(FREE_SURF==0){
	        for (i=1;i<=NX;i++){
		    for (j=1;j<=NPML;j++){
		        A[j][i] *= absorb_coeff[j][i];
			if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}
		    }
	        }
            }

	     /* bottom boundary */
	     for (i=1;i<=NX;i++){
		    for (j = FSSHIFT + NY0 + 1;j <= NY;j++){
			    A[j][i] *= absorb_coeff[j][i];
			    if(isnan(A[j][i])==1){A[j][i]=A[j-1][i];}
		    }
	     }

	    sprintf(modfile,"grad_damp_PML.bin");
	    writemod(modfile,absorb_coeff,3); 

	    free_vector(coeff,1,NPML);	
	    free_matrix(absorb_coeff,1,NY,1,NX);
	
	}
			
}
