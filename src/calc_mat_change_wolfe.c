/*------------------------------------------------------------------------
 *   update material parameter for Wolfe line search 
 *   
 *   Daniel Koehn
 *   last update 14.12.2015
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void calc_mat_change_wolfe(float  **  Hgrad, float **  vp, float **  vp_old, float eps_scale, int itest){

	/* global variables */
	extern int NX, NY, MYID;
	extern char INV_MODELFILE[STRING_SIZE];
	extern float MAT1_LOW, MAT1_UP;

	/* local variables */
	int i, j;
	char modfile[STRING_SIZE];
        float maxgrad, maxvp;

        /* constant step length for debugging */
        /*maxgrad = maximum_m(Hgrad,NX,NY);
        maxvp = maximum_m(Vp,NX,NY);
        eps_scale = 0.01 * maxvp/maxgrad;*/	

	/* loop over local grid */
	for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){		
		    
		    /* P-wave velocity */
		    vp[j][i] = vp_old[j][i] + eps_scale * Hgrad[j][i]; 	
		  
		    /* apply hard constraints */
	      	    if(vp[j][i]<MAT1_LOW){
	               vp[j][i] = vp_old[j][i];
	            }
		      
		    if(vp[j][i]>MAT1_UP){
		       vp[j][i] = vp_old[j][i];
		    }
		      
		      
		    /* P-wave velocity should not be smaller than zero */
		    if(vp[j][i]<0.0){
		       vp[j][i] = vp_old[j][i];
		    }
  
		    /*if(itest==0){
		       Vp[j][i] = Vpnp1[j][i]; 
	            } */
		                
		}
	}
	
	if((itest==0)&&(MYID==0)){
	   sprintf(modfile,"%s_vp.bin",INV_MODELFILE);
	   writemod_true(modfile,vp,3);
	}
                                                
}

