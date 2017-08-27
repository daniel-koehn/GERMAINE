/*------------------------------------------------------------------------
 *   update all material parameter classes for Wolfe line search 
 *   
 *   Daniel Koehn
 *   last update 22.08.2017
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void calc_mat_change_wolfe_multi_para(float  **  Hgrad, float **  vp, float **  vp_old, float eps_scale, int para_index){

	/* global variables */
	extern int NX, NY, MYID;
	extern char INV_MODELFILE[STRING_SIZE];
	extern float MAT1_LOW, MAT1_UP, MAT2_LOW, MAT2_UP;

	/* local variables */
	int i, j;
	char modfile[STRING_SIZE];
        float maxgrad, maxvp, vp_low, vp_up;

	if(para_index==1){
	    vp_low = MAT1_LOW;
	    vp_up = MAT1_UP;
	}

	if(para_index==2){
	    vp_low = MAT2_LOW;
	    vp_up = MAT2_UP;
	}


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
	      	    if(vp[j][i] < vp_low){
	               vp[j][i] = vp_old[j][i];
	            }
		      
		    if(vp[j][i] > vp_up){
		       vp[j][i] = vp_old[j][i];
		    }
		      		      
		    /* P-wave velocity should not be smaller than zero */
		    if(vp[j][i]<0.0){
		       vp[j][i] = vp_old[j][i];
		    }
		                
		}
	}	
                                                
}

