/*
 * Assemble gradients from all shots  
 *
 * Daniel Koehn
 * Kiel, 20.08.2017
 *
 */

#include "fd.h"

void ass_grad_TE(struct fwiTE *fwiTE, struct waveAC *waveAC, struct matTE *matTE, float ** grad_shot_sigma, float ** grad_shot_epsilon, float **srcpos,  int nshots, int **recpos, int ntr, int ishot){

        /* global variables */
	extern int NX, NY, SWS_TAPER_CIRCULAR_PER_SHOT;
	extern float S;

	/* local variables */
	int i, j;
	complex float wien, Omega, Omega2;

	wien = 	(*waveAC).stfr + (*waveAC).stfi * I;	
	Omega = 2.0 * M_PI * (*waveAC).freq - (I * S);
	Omega2 = cpowf(((2.0*M_PI*(*waveAC).freq) - (I * S)),2.0);

	/* printf("abs(wien) = %e + i %e \n",creal(cabsf(wien)),cimag(cabsf(wien))); */
	
	/* assemble sigma gradient for one shot */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

               grad_shot_sigma[j][i] =  creal( (wien/cabsf(wien)) * ((*fwiTE).forwardr[j][i] + (*fwiTE).forwardi[j][i] * I)  
				       * (Omega * I) * ((*waveAC).pr[j][i] + (*waveAC).pi[j][i] * I) );

	       grad_shot_epsilon[j][i] =  creal(Omega2 * (wien/cabsf(wien)) * ((*fwiTE).forwardr[j][i] + (*fwiTE).forwardi[j][i] * I)  
				       * ((*waveAC).pr[j][i] + (*waveAC).pi[j][i] * I) );		    

	    }
	}

	/* apply taper at source positions */
	if(SWS_TAPER_CIRCULAR_PER_SHOT==1){
	  taper_grad_shot(grad_shot_sigma,srcpos,nshots,recpos,ntr,ishot);
	  taper_grad_shot(grad_shot_epsilon,srcpos,nshots,recpos,ntr,ishot);
	}

	if(SWS_TAPER_CIRCULAR_PER_SHOT==2){
	  taper_grad_shot1(grad_shot_sigma,srcpos,nshots,recpos,ntr,ishot);
	  taper_grad_shot1(grad_shot_epsilon,srcpos,nshots,recpos,ntr,ishot);
	}
	
	/* assemble gradients from all shots */
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){

	       (*fwiTE).grad_sigma[j][i] += grad_shot_sigma[j][i];
	       (*fwiTE).grad_epsilon[j][i] += grad_shot_epsilon[j][i];
		    
	    }
	}
	

}



