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
	extern float S;

	/* local variables */
	int i, j;
	complex float wien, Omega2;

	wien = 	(*waveAC).stfr + (*waveAC).stfi * I;	
	Omega2 = cpowf(((2.0*M_PI*(*waveAC).freq) - (I * S)),2.0);

	/* printf("abs(wien) = %e + i %e \n",creal(cabsf(wien)),cimag(cabsf(wien))); */
	
	/* assemble gradient for one shot */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

               grad_shot[j][i] = - 2.0 * creal(Omega2 * (wien/cabsf(wien)) * ((*fwiAC).forwardr[j][i] + (*fwiAC).forwardi[j][i] * I)  
				       * (1.0 / ((*matAC).rho[j][i] * (*matAC).vp[j][i] * (*matAC).vp[j][i] * (*matAC).vp[j][i])) * ((*waveAC).pr[j][i] + (*waveAC).pi[j][i] * I) );
		    
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



