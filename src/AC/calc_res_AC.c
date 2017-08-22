/*------------------------------------------------------------------------
 *  Calculate FD residuals and objective function                           
 *  
 *  D. Koehn
 *  Kiel, 23/06/2016
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float calc_res_AC(struct waveAC *waveAC, int ntr, int ishot, int nstage, int nfreq){

        /* global variables */
	extern int STF_INV, INVMAT, MISFIT, NSHOTS;
	extern char DATA_DIR[STRING_SIZE];

	/* local variables */
	int i, index;
        float l2;
        complex float res, wiennom, wiendenom, wien;

	if(STF_INV==1){           

	   l2 = 0.0;

           /* initiate wiener deconvolution */
           wiennom = 0.0 + 0.0 * I;
           wiendenom = 0.0 + 0.0 * I;

	   for(i=1;i<=ntr;i++){  

	       index = i + ntr * (ishot-1) + ntr * NSHOTS * (nfreq-1);
		
	       /* estimate nominator and denominator of Wiener deconvolution */
	    
	       wiennom += conj((*waveAC).precr[index] + (*waveAC).preci[index] * I) * ((*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I);
	       wiendenom += conj((*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I) * ((*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I);
	          	    
	   }

           /* estimate STF */        
           wien = wiennom / wiendenom;	

	}

	if(STF_INV==0){
	    wien = 1.0 + 0.0 * I;
        }

        (*waveAC).stfr = creal(wien);
        (*waveAC).stfi = cimag(wien);

	/* printf("stfr = %e \t stfi = %e \n",(*waveAC).stfr,(*waveAC).stfi); */

	for(i=1;i<=ntr;i++){

	    index = i + ntr * (ishot-1) + ntr * NSHOTS * (nfreq-1);

	    /* calculate complex data residuals ... */
    	    /* ... for FDFD-FWI */
	    if(INVMAT==1){
		
		if(MISFIT==1){ /* L2-norm */
	            res = (((*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I) - (wien*((*waveAC).precr[index] + (*waveAC).preci[index] * I)))/cabsf(wien);
		}

		if(MISFIT==2){ /* logarithmic L2-norm (Shin and Min 2006) */
	            res = clogf((*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I) - clogf((*waveAC).precr[index] + (*waveAC).preci[index] * I);
		}

		if(MISFIT==3){ /* phase-only residuals (Bednar et al. 2007) */
	            res = cargf((*waveAC).precr[index] + (*waveAC).preci[index] * I) - cargf((*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I);
		}

	    }

    	    /* ... for FDFD-RTM */
	    if(INVMAT==2){res = (*waveAC).pobsr[index] + (*waveAC).pobsi[index] * I;}

	        (*waveAC).presr[i] = creal(res);
	        (*waveAC).presi[i] = cimag(res);

	    if(((*waveAC).pobsr[index]<1e-20)&&((*waveAC).pobsi[index]<1e-20)){
                (*waveAC).presr[i] = 0.0;
	    	(*waveAC).presi[i] = 0.0;
	    }

            /* calculate objective function */
	    l2 += creal(res * conj(res));
 
	}

	/*printf("l2 = %e \n",l2);*/

        return l2;
}
