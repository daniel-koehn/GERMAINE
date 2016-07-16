/*
 * Assemble gradient from all shots  
 *
 * Daniel Koehn
 * Kiel, 29.06.2016
 *
 */

#include "fd.h"

void ass_grad_AC(struct fwiAC *fwiAC, struct waveAC *waveAC, struct matAC *matAC, float ** grad_shot, float **srcpos,  int nshots, int **recpos, int ntr, int ishot){

        /* global variables */
	extern int NX, NY, SWS_TAPER_CIRCULAR_PER_SHOT;

	/* local variables */
	int i, j;
	complex float wien;

	wien = 	(*fwiAC).stfr + (*fwiAC).stfi * I;	

	// printf("abs(wien) = %e + i %e \n",creal(fabs(wien)),cimag(fabs(wien)));
	
	/* assemble gradient for one shot */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

               grad_shot[j][i] = - 2.0 * (*waveAC).omega2 * creal( wien/cabsf(wien) * ((*fwiAC).forwardr[j][i] + (*fwiAC).forwardi[j][i] * I)  
				       * (1.0 / ((*matAC).vp[j][i] * (*matAC).vp[j][i] * (*matAC).vp[j][i])) * ((*waveAC).pr[j][i] + (*waveAC).pi[j][i] * I) );
		    
	    }
	}

	/* apply taper at source positions */
	if(SWS_TAPER_CIRCULAR_PER_SHOT==1){
	  taper_grad_shot(grad_shot,srcpos,nshots,recpos,ntr,ishot);
	}

	if(SWS_TAPER_CIRCULAR_PER_SHOT==2){
	  taper_grad_shot1(grad_shot,srcpos,nshots,recpos,ntr,ishot);
	}
	
	/* assemble gradient from all shots */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

	       (*fwiAC).grad[j][i] += grad_shot[j][i];
		    
	    }
	}
	

}



