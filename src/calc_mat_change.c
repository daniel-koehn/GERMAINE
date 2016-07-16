/*------------------------------------------------------------------------
 *   calculate test step length for material parameter update
 *   
 *   Daniel Koehn
 *   last update 11.12.2015
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
float calc_mat_change(float  **  waveconv, float **  pi, float **  pinp1, int iter, float eps_scale, int itest, int nfstart){

	/*--------------------------------------------------------------------------*/
	FILE *FP1;
	/* extern variables */
	extern float DH, DT;
	extern float EPSILON;
	extern int NX, NY, MYID;
	
	extern char INV_MODELFILE[STRING_SIZE];
	
	extern float VPUPPERLIM, VPLOWERLIM;

	/* local variables */

	float Zp, eps_true, pimax, gradmax, epsilon1;
	int i, j;
	char modfile[STRING_SIZE];
	
	/* invert for Zp and Zs */
	/* ------------------------------------------------------------------------------------ */	
	
	/* find maximum of Zp and gradient waveconv */	

	pimax = 0.0;
	gradmax = 0.0;
	
	    for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
		
		Zp = pi[j][i];
		
		
		if(Zp>pimax){pimax=Zp;}
		
		if((i*j == 1) || (gradmax == 0.0)) {
				gradmax = fabs(waveconv[j][i]);		
		} else {		
		   if(fabs(waveconv[j][i]) > gradmax){
		      gradmax = fabs(waveconv[j][i]);
		   }
		                               
		}
	}}	
	
	
	/* calculate scaling factor for the gradient */
        /* ----------------------------------------- */
/* 	EPSILON = eps_scale * (pimax/gradmax); */
	EPSILON = 0.005 * (pimax/gradmax);
		
        if(MYID==0){
	  printf("|Vpmax| = %e \t |gradmax| = %e \n",pimax,gradmax);
	}	
      
      /* save true step length */
      eps_true = EPSILON;
	
	/* loop over local grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){		
		    
		    /* P-wave velocity */
		    pinp1[j][i] = pi[j][i] - EPSILON*waveconv[j][i]; 	
		  
		    /* apply hard constraints */
	      	    if(pinp1[j][i]<VPLOWERLIM){
	               pinp1[j][i] = pi[j][i];
	            }
		      
		    if(pinp1[j][i]>VPUPPERLIM){
		       pinp1[j][i] = pi[j][i];
		    }
		      
		      
		    /* P-wave velocity should not be smaller than zero */
		    if(pinp1[j][i]<0.0){
		       pinp1[j][i] = pi[j][i];
		    }
  
		    if(itest==0){
		       pi[j][i] = pinp1[j][i]; 
	            } 
		                
		}
	}
	
	if(itest==0){
	   sprintf(modfile,"%s_vp.bin",INV_MODELFILE);
	   writemod(modfile,pinp1,3);
	}
        
        if((itest==0)&&(iter==nfstart)){
	   sprintf(modfile,"%s_vp_it_%d.bin",INV_MODELFILE,iter);
	   writemod(modfile,pinp1,3);                                        
	}
                                                

	return eps_true;
}

