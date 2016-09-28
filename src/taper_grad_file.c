/*------------------------------------------------------------------------
 * Apply taper function from file  
 * 
 * Daniel Koehn
 * Kiel, 03.07.2016
 * ----------------------------------------------------------------------*/

#include "fd.h"

void taper_grad_file(float ** grad){

        /* global variables */
	extern int NX, NY, NX0, NY0, NPML, FREE_SURF, FSSHIFT;
	extern int SWS_TAPER_FILE;
	
        /* local variables */
	int i, j;
	float grad_tap;
	
	if (SWS_TAPER_FILE){	

	   FILE *fp_taper;	   

	   fp_taper=fopen("taper.bin","r");
	
	   for (i=1;i<=NX0;i++){
	        for (j=1;j<=NY0;j++){

		    fread(&grad_tap, sizeof(float), 1, fp_taper); 
		    grad[j+FSSHIFT][i+NPML] *= grad_tap;		

		}
	   }
	   
	   /* extend model */
	   
	   /* left boundary */
	   for (i=1;i<=NPML;i++){
		   for (j=1;j<=NY;j++){
			   grad[j][i] = grad[j][NPML+1];				
		   }
	   }
	
	   /* right boundary */
	   for (i = NPML + NX0 + 1;i <= NX;i++){
		   for (j=1;j<=NY;j++){
			   grad[j][i] = grad[j][NPML + NX0];				
		   }
	   }

	   /* top boundary */
	   if(FREE_SURF==1){
	       for (i=1;i<=NX;i++){
	           for (j=1;j<=NPML;j++){
		       grad[j][i] = grad[NPML+1][i];				
		   }
	       }
	   }

	   /* bottom boundary */
	   for (i=1;i<=NX;i++){
		   for (j = FSSHIFT + NY0 + 1;j <= NY;j++){
			grad[j][i] = grad[FSSHIFT+NY0][i];				
		   }
	   } 
		 
	   
	   fclose(fp_taper);   

	}

}
